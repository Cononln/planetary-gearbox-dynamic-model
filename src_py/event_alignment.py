"""Event-constrained alignment pilot. All time shifts are in seconds externally.

The relative path uses a zero-channel-mean gauge. Common residual shifts do
not uniquely identify speed error. No encoder is passed to event inference.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
from scipy import signal, ndimage, optimize, io


def interval(t, span):
    return (t >= span[0]) & (t < span[1])


def rational_resample(x, source_fs, target_fs):
    ratio = Fraction(int(target_fs), int(source_fs))
    return signal.resample_poly(x, ratio.numerator, ratio.denominator, axis=0).astype(np.float32)


def pulse_edges(x, fs):
    lo, hi = np.percentile(x[::101], [2, 98])
    if hi - lo <= 1e-9:
        return np.empty(0)
    threshold = (lo + hi) / 2
    ix = np.flatnonzero((x[1:] >= threshold) & (x[:-1] < threshold)) + 1
    if not len(ix):
        return np.empty(0)
    ix = ix[np.r_[True, np.diff(ix) > int(fs * .00015)]]
    return ix / fs


def load_record(path, cfg):
    a = io.loadmat(path, variable_names=['Data', 'SampleRate'])
    x = a['Data']
    fs = float(np.asarray(a['SampleRate']).squeeze())
    if x.ndim != 2 or x.shape[1] != 7 or fs <= 0:
        raise ValueError('Unexpected data schema')
    if not np.isfinite(x).all():
        raise ValueError('Non-finite samples must be audited before inference')
    enc = [pulse_edges(x[:, j], fs) for j in cfg['encoder_evaluation_columns']]
    high = rational_resample(x[:, cfg['vibration_columns']], fs, cfg['working_rates_hz']['impulse'])
    del a, x
    fhi = cfg['working_rates_hz']['impulse']
    flo = cfg['working_rates_hz']['motion']
    low = rational_resample(high, fhi, flo)
    tc = np.arange(len(low)) / flo
    mean = low[interval(tc, cfg['split_seconds']['calibration'])].mean(0)
    high -= mean
    low -= mean
    return high, low, fs, enc


def estimate_motion(low, cfg):
    fs = cfg['working_rates_hz']['motion']
    t = np.arange(len(low)) / fs
    cal = interval(t, cfg['split_seconds']['calibration'])
    y = low[cal]
    f = np.fft.rfftfreq(len(y), 1 / fs)
    a = np.sqrt(np.mean(np.abs(np.fft.rfft(y * np.hanning(len(y))[:, None], axis=0)) ** 2, axis=1))
    limits = cfg['motion']['mesh_search_hz']
    use = np.flatnonzero((f >= limits[0]) & (f <= limits[1]))
    k = use[np.argmax(a[use])]
    logs = np.log(np.maximum(a[k-1:k+2], 1e-30))
    df = .5 * (logs[0]-logs[2]) / (logs[0]-2*logs[1]+logs[2])
    fm = f[k] + np.clip(df, -1, 1) * fs / len(y)
    num = np.zeros(len(t)-1)
    den = np.zeros(len(t)-1)
    for h in cfg['motion']['harmonics']:
        bw = 30 + 4 * (h-1)
        taps = signal.firwin(2*round(.25*fs)+1, bw, fs=fs, window=('kaiser', 6))
        osc = np.exp(2j*np.pi*h*fm*t)
        for j in range(3):
            z = 2 * signal.fftconvolve(low[:, j] * osc.conj(), taps, mode='same') * osc
            amp = np.abs(z)
            r = amp / max(np.median(amp[cal]), 1e-20)
            w = r*r/(1+r*r)
            weight = np.minimum(w[1:], w[:-1])
            inc = np.angle(z[1:]*z[:-1].conj()*np.exp(-2j*np.pi*h*fm/fs)) / h
            num += inc * weight
            den += weight
    inc = ndimage.uniform_filter1d(num / np.maximum(den, 1e-20), 2*int(.0015*fs//2)+1, mode='nearest')
    drift = ndimage.uniform_filter1d(np.r_[0, np.cumsum(inc)], 2*int(.012*fs//2)+1, mode='nearest')
    slope = np.polyfit(t[cal], drift[cal], 1)[0]
    phase = 2*np.pi*fm*t + drift - slope*t
    if np.any(np.diff(phase) <= 0):
        raise ValueError('Nonmonotonic vibration-derived motion phase')
    return dict(time=t, phase=phase, fm=float(fm), slope_bias_hz=float(slope/(2*np.pi)))


def make_evidence(high, motion, band, cfg):
    fhi = cfg['working_rates_hz']['impulse']
    flo = cfg['working_rates_hz']['motion']
    length = cfg['transient']['analytic_fir_taps']
    n = np.arange(length) - (length-1)/2
    fc, half = np.mean(band), (band[1]-band[0])/2
    h = 2*signal.firwin(length, half, fs=fhi, window=('kaiser', 8))*np.exp(2j*np.pi*fc*n/fhi)
    z = np.empty(high.shape, np.complex64)
    energy = np.empty((len(motion['time']), 3), np.float32)
    for j in range(3):
        z[:, j] = signal.fftconvolve(high[:, j], h, mode='same')
        e = ndimage.uniform_filter1d(np.abs(z[:, j])**2, max(1, round(.0005*fhi)))
        energy[:, j] = rational_resample(e, fhi, flo)[:len(energy)]
    return z.real.copy(), contrast_from_energy(energy, motion, cfg), energy


def contrast_from_energy(energy, motion, cfg, events=None):
    phase = motion['phase']
    t = motion['time']
    cal = interval(t, cfg['split_seconds']['calibration'])
    d = np.zeros_like(energy)
    for j in range(3):
        neighbors = []
        for q in cfg['transient']['mesh_neighbor_offsets']:
            values = np.interp(phase+2*np.pi*q, phase, energy[:, j], left=np.nan, right=np.nan)
            if events is not None:
                query_t = np.interp(phase+2*np.pi*q, phase, t, left=np.nan, right=np.nan)
                ix = np.searchsorted(events['t'], query_t+.012)-1
                near = (ix >= 0) & (events['t'][np.clip(ix, 0, len(events['t'])-1)] >= query_t-.024)
                values[near] = np.nan
            neighbors.append(values)
        neighbors = np.stack(neighbors)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            center = np.nanmedian(neighbors, axis=0)
            mad = 1.4826*np.nanmedian(np.abs(neighbors-center), axis=0)
        floor = cfg['transient']['contrast_denominator_floor_fraction'] * np.median(energy[cal, j])
        d[:, j] = np.nan_to_num((energy[:, j]-center)/np.maximum(mad+floor, 1e-20))
        d[np.isfinite(neighbors).sum(axis=0) < 4, j] = 0
    return d


def event_times(motion, eta, start, end, zp=31, order_factor=1.0):
    phase = motion['phase'] / zp * order_factor
    p0, p1 = np.interp([start, end], motion['time'], phase)
    n = np.arange(np.ceil((p0-eta)/(2*np.pi)), np.floor((p1-eta)/(2*np.pi))+1)
    times = np.interp(eta + 2*np.pi*n, phase, motion['time'])
    return times, n.astype(int)


def schedule(d, motion, cfg, order_factor=1.0):
    cal = cfg['split_seconds']['calibration']
    t = motion['time']
    def score(eta):
        tt, _ = event_times(motion, eta, *cal, order_factor=order_factor)
        yy = np.stack([np.interp(tt, t, d[:, j]) for j in range(3)], axis=1)
        support = np.sort(yy, axis=1)[:, 1]
        return float(np.mean(np.clip(support, -3, 10)))
    grid = np.deg2rad(np.arange(0, 360, cfg['events']['initial_phase_grid_deg']))
    values = np.array([score(eta) for eta in grid])
    eta_a = grid[np.argmax(values)]
    offsets = np.arange(.35, .65001, .01)*2*np.pi
    eta_b = eta_a + offsets[np.argmax([score(eta_a+offset) for offset in offsets])]
    rows = []
    for c, eta in enumerate([eta_a, eta_b]):
        tt, number = event_times(motion, eta, 2.5, 59.3, order_factor=order_factor)
        for time, index in zip(tt, number):
            psi = np.interp(time, t, motion['phase'])/cfg['gear_teeth']['ring']
            rows.append((time, psi, c, index))
    rows.sort()
    arr = np.asarray(rows)
    return dict(t=arr[:, 0], psi=arr[:, 1], branch=arr[:, 2].astype(int), recurrence=arr[:, 3].astype(int),
                origin_rad=np.array([eta_a, eta_b]), origin_scores=np.array([score(eta_a), score(eta_b)]),
                order_factor=order_factor)


def candidates(events, d, motion, cfg):
    fs = cfg['working_rates_hz']['motion']
    maxk = cfg['events']['max_candidates_per_channel']
    bound = cfg['events']['search_half_width_s']
    threshold = cfg['transient']['initial_contrast_threshold']
    offsets = np.full((len(events['t']), 3, maxk), np.nan)
    scores = np.full_like(offsets, -np.inf)
    for n, time in enumerate(events['t']):
        a, b = int(np.ceil((time-bound)*fs)), int(np.floor((time+bound)*fs))+1
        for j in range(3):
            peaks, _ = signal.find_peaks(d[a:b, j], height=threshold,
                                        distance=max(1, int(cfg['events']['peak_min_separation_s']*fs)))
            peaks = peaks[np.argsort(d[a+peaks, j])[-maxk:][::-1]]
            for k, p in enumerate(peaks):
                index = a+p
                y = d[index-1:index+2, j].astype(float)
                sub = .5*(y[0]-y[2]) / min(y[0]-2*y[1]+y[2], -1e-12)
                offsets[n, j, k] = (index+np.clip(sub, -.5, .5))/fs-time
                scores[n, j, k] = d[index, j]
    return dict(offsets=offsets, scores=scores)


def path_design(psi, branch, variant='full'):
    psi, branch = np.asarray(psi), np.asarray(branch)
    basis = np.zeros((len(psi), 3, 20))
    for n, (angle, c) in enumerate(zip(psi, branch)):
        start = 10*int(c)
        basis[n, :, start:start+2] = [[1, 0], [0, 1], [-1, -1]]
        alpha = np.arange(3)*2*np.pi/3
        for k in [1, 2]:
            co, si = np.cos(k*(angle-alpha)), np.sin(k*(angle-alpha))
            basis[n, :, start+2*k:start+2*k+2] = np.stack([co-co.mean(), si-si.mean()], axis=1)
        u = np.array([np.cos(angle), np.sin(angle)])
        basis[n, 0, start+6:start+8] = u
        basis[n, 1, start+8:start+10] = u
        basis[n, 2, start+6:start+8] = -u
        basis[n, 2, start+8:start+10] = -u
    if variant == 'independent':
        basis[:, :, [2, 3, 4, 5, 12, 13, 14, 15]] = 0
    return basis


def candidate_selection(cand, predicted_ms, excluded=None):
    offsets = cand['offsets']*1000
    score = cand['scores']
    weights = np.exp(np.minimum(score, 10)/2)
    norm = np.exp(1.5)+weights.sum(axis=2)
    p = weights / norm[:, :, None]
    missing = -np.log(np.exp(1.5)/norm)
    residual = offsets-predicted_ms[:, :, None]
    a = np.abs(residual)
    rho = np.where(a <= 1.5, .5*residual**2, 1.5*(a-.75))
    cost = np.where(np.isfinite(offsets), rho-np.log(np.maximum(p, 1e-30)), np.inf)
    best = np.argmin(cost, axis=2)
    value = np.take_along_axis(offsets, best[:, :, None], axis=2)[:, :, 0]
    quality = np.take_along_axis(cost, best[:, :, None], axis=2)[:, :, 0]
    valid = quality < missing
    if excluded is not None:
        valid[:, excluded] = False
    value[~valid] = np.nan
    return value, valid


def smooth_matrix(t):
    n = len(t)
    if n < 3:
        return np.zeros((0, n))
    dt = np.diff(t)
    ref = np.median(dt)
    result = np.zeros((n-2, n))
    for k in range(n-2):
        result[k, k:k+3] = [ref/dt[k], -ref/dt[k]-ref/dt[k+1], ref/dt[k+1]]
    return result


def fit_delay_model(events, cand, span, smooth=1.0, variant='full'):
    use = interval(events['t'], span)
    t = events['t'][use]
    basis = path_design(events['psi'][use], events['branch'][use], variant)
    small = {k: v[use] for k, v in cand.items()}
    raw = small['offsets'][:, :, 0]*1000
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        common = np.nanmedian(raw, axis=1)
    common = np.clip(np.nan_to_num(common), -4, 4)
    coeff = np.zeros(20)
    n = len(t)
    sm = smooth_matrix(t)
    grid = path_design(np.tile(np.linspace(0, 2*np.pi, 73), 2), np.repeat([0, 1], 73), variant).reshape(-1, 20)
    matrix = np.c_[np.zeros((len(grid), n)), grid]
    linear = optimize.LinearConstraint(matrix, -4, 4)
    penalty = np.tile([1e-5, 1e-5, .03, .03, .48, .48, .3, .3, .3, .3], 2)
    if variant == 'independent':
        penalty[[6, 7, 8, 9, 16, 17, 18, 19]] = .03
    if variant == 'no_path_smoothing':
        penalty[:] = 1e-6
    history = []
    previous_assignment = None
    for iteration in range(10):
        prediction = common[:, None]+basis@coeff
        values, valid = candidate_selection(small, prediction)
        event_valid = valid.sum(axis=1) >= 2
        valid &= event_valid[:, None]
        ni, ji = np.where(valid)
        if len(ni) < 12:
            return dict(coeff=np.zeros(20), smooth=smooth, variant=variant, history=history,
                        success=False, reason='insufficient_supported_events', calibration_coverage=float(event_valid.mean()))
        xx, yy = basis[ni, ji], values[ni, ji]
        def objective(z):
            d, beta = z[:n], z[n:]
            r = d[ni]+xx@beta-yy
            ar = np.abs(r)
            rho = np.where(ar < 1.5, .5*r*r, 1.5*(ar-.75))
            loss = rho.sum()/(3*n)
            clip = np.clip(r, -1.5, 1.5)/(3*n)
            gd = np.bincount(ni, weights=clip, minlength=n)
            gb = xx.T@clip
            if len(sm):
                sd = sm@d
                loss += smooth*np.mean(sd*sd)
                gd += 2*smooth*(sm.T@sd)/len(sm)
            loss += np.sum(penalty*beta*beta)/20 + 1e-7*np.sum(d*d)
            gb += 2*penalty*beta/20
            gd += 2e-7*d
            return loss, np.r_[gd, gb]
        z0 = np.r_[common, coeff]
        result = optimize.minimize(objective, z0, jac=True, method='SLSQP',
                                   bounds=[(-4, 4)]*n+[(-10, 10)]*20, constraints=[linear],
                                   options=dict(maxiter=1000, ftol=1e-7))
        if not result.success or np.max(np.abs(grid@result.x[n:])) > 4.001:
            return dict(coeff=coeff, smooth=smooth, variant=variant, history=history, success=False,
                        reason=str(result.message), solver_iterations=int(result.nit),
                        solver_loss=float(result.fun), calibration_coverage=float(event_valid.mean()))
        common, coeff = result.x[:n], result.x[n:]
        history.append(dict(iteration=iteration, loss=float(result.fun), accepted=int(event_valid.sum()),
                            solver_iterations=int(result.nit), parameter_update_ms=float(np.max(np.abs(result.x-z0)))))
        if previous_assignment is not None and np.array_equal(np.nan_to_num(values), previous_assignment) and history[-1]['parameter_update_ms'] < .01:
            break
        previous_assignment = np.nan_to_num(values)
    covered = np.unique(np.floor(np.mod(events['psi'][use][event_valid], 2*np.pi)/(np.pi/3)).astype(int))
    okay = len(covered) >= 4
    return dict(coeff=coeff, smooth=smooth, variant=variant, history=history, success=okay,
                reason='okay' if okay else 'insufficient_carrier_coverage', calibration_coverage=float(event_valid.mean()))


def predict_delays(events, cand, model, spans, excluded=None):
    path_ms = path_design(events['psi'], events['branch'], model['variant'])@model['coeff']
    shift = np.zeros_like(path_ms)
    accepted = np.zeros(len(path_ms), bool)
    measured = np.full_like(path_ms, np.nan)
    common_all = np.zeros(len(path_ms))
    for span in spans:
        use = interval(events['t'], span)
        ids = np.flatnonzero(use)
        small = {k: v[use] for k, v in cand.items()}
        raw = small['offsets'][:, :, 0]*1000-path_ms[use]
        if excluded is not None:
            raw[:, excluded] = np.nan
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            common = np.clip(np.nan_to_num(np.nanmedian(raw, axis=1)), -4, 4)
        sm = smooth_matrix(events['t'][use])
        for _ in range(6):
            values, valid = candidate_selection(small, common[:, None]+path_ms[use], excluded)
            ev = valid.sum(axis=1) >= 2
            valid &= ev[:, None]
            count = valid.sum(axis=1).astype(float)
            target = np.nansum(np.where(valid, values-path_ms[use], np.nan), axis=1)
            matrix = np.diag(count+1e-5) + model['smooth']*(sm.T@sm)
            common = np.clip(np.linalg.solve(matrix, target), -4, 4)
        shift[ids] = (common[:, None]+path_ms[use])/1000
        accepted[ids] = ev
        measured[ids] = values/1000
        common_all[ids] = common/1000
    if not model['success']:
        accepted[:] = False
        shift[:] = 0
    shift[~accepted] = 0
    return dict(shift=shift, accepted=accepted, measured=measured, path=path_ms/1000, common=common_all)


def sample_sinc(x, query_samples, half=16):
    """Windowed-sinc interpolation; constant integer shifts are exact."""
    q = np.asarray(query_samples)
    integer = np.floor(q).astype(np.int64)
    result = np.zeros(q.shape)
    weight_sum = np.zeros(q.shape)
    for k in range(-half+1, half+1):
        ix = integer+k
        distance = q-ix
        weight = np.sinc(distance)*np.sinc(distance/half)
        okay = (ix >= 0) & (ix < len(x))
        result += weight*np.where(okay, x[np.clip(ix, 0, len(x)-1)], 0)
        weight_sum += weight
    return result/np.maximum(weight_sum, 1e-12)


def extract_events(band, fs, events, shifts, core=(-.004, .016)):
    u = np.arange(round(core[0]*fs), round(core[1]*fs))/fs
    out = np.empty((len(events['t']), 3, len(u)), np.float32)
    for j in range(3):
        query = (events['t'][:, None]+shifts[:, j, None]+u)*fs
        out[:, j] = sample_sinc(band[:, j], query)
    return out, u


def make_templates(w, events, cal):
    templates = {}
    bins = (np.mod(events['psi'], 2*np.pi)/(np.pi/3)).astype(int)
    for c in [0, 1]:
        for b in range(6):
            use = interval(events['t'], cal) & (events['branch'] == c) & (bins == b)
            if use.sum() < 3:
                continue
            for j in range(3):
                template = np.median(w[use, j], axis=0)
                templates[(c, b, j)] = template-template.mean()
    return templates


def local_xcorr(band, fs, events, templates, accepted):
    shifts = np.zeros((len(events['t']), 3))
    bins = (np.mod(events['psi'], 2*np.pi)/(np.pi/3)).astype(int)
    width, search = round(.020*fs), round(.008*fs)
    for n, (t0, c, b) in enumerate(zip(events['t'], events['branch'], bins)):
        if not accepted[n]:
            continue
        center = round((t0-.004)*fs)
        for j in range(3):
            template = templates.get((int(c), int(b), j))
            if template is None or np.linalg.norm(template) < 1e-15:
                continue
            source = band[center-search:center+width+search, j].astype(float)
            score = signal.correlate(source, template, mode='valid', method='fft')
            sums = np.r_[0, np.cumsum(source)]
            squares = np.r_[0, np.cumsum(source*source)]
            variance = squares[width:]-squares[:-width]-(sums[width:]-sums[:-width])**2/width
            score /= np.maximum(np.sqrt(np.maximum(variance, 0))*np.linalg.norm(template), 1e-20)
            k = int(np.argmax(score))
            extra = 0.0
            if 0 < k < len(score)-1:
                denominator = score[k-1]-2*score[k]+score[k+1]
                extra = np.clip(.5*(score[k-1]-score[k+1])/min(denominator, -1e-15), -.5, .5)
            shifts[n, j] = (k-search+extra)/fs + center/fs-(t0-.004)
    return np.clip(shifts, -.008, .008)


def reconstruct_block(band, fs, events, shifts, accepted, span):
    a, b = [round(s*fs) for s in span]
    t = np.arange(a, b)/fs
    displacement = np.zeros((b-a, 3))
    dt = np.diff(events['t'])
    gaps = np.minimum(np.r_[np.inf, dt], np.r_[dt, np.inf])
    rejected = 0
    for n in np.flatnonzero(accepted & (events['t'] > span[0]-.1) & (events['t'] < span[1]+.1)):
        length = min(.020, (.8*gaps[n]-.020)/2)
        if length <= 0 or np.max(np.abs(shifts[n]))*np.pi/(2*length) > .8:
            rejected += 1
            continue
        rel = t-events['t'][n]
        gate = np.zeros(len(t))
        core = (rel >= -.004) & (rel < .016)
        left = (rel >= -.004-length) & (rel < -.004)
        right = (rel >= .016) & (rel < .016+length)
        gate[core] = 1
        gate[left] = .5*(1-np.cos(np.pi*(rel[left]+.004+length)/length))
        gate[right] = .5*(1+np.cos(np.pi*(rel[right]-.016)/length))
        displacement += gate[:, None]*shifts[n]
    derivative = 1+np.diff(displacement, axis=0)*fs
    if derivative.min(initial=1) < .199:
        raise ValueError('Time map failed monotonicity check')
    output = np.empty((len(t), 3), np.float32)
    for j in range(3):
        output[:, j] = sample_sinc(band[:, j], np.arange(a, b)+displacement[:, j]*fs)
    return output, float(derivative.min(initial=1)), rejected


def cross_block_metrics(w, events, mask, blocks):
    bins = (np.mod(events['psi'], 2*np.pi)/(np.pi/3)).astype(int)
    rows = []
    for c in [0, 1]:
        for b in range(6):
            for j in range(3):
                means, raw_means, energies = [], [], []
                for span in blocks:
                    use = mask & interval(events['t'], span) & (events['branch'] == c) & (bins == b)
                    yy = w[use, j].astype(float)
                    if len(yy) < 2:
                        means.append(None); raw_means.append(None); energies.append(np.nan)
                        continue
                    yy -= yy.mean(axis=1, keepdims=True)
                    yy *= signal.windows.tukey(yy.shape[1], .2)[None, :]
                    norms = np.linalg.norm(yy, axis=1)
                    if np.any(norms <= 1e-15):
                        means.append(None); raw_means.append(None); energies.append(np.nan)
                        continue
                    means.append(np.mean(yy/norms[:, None], axis=0))
                    raw_means.append(np.mean(yy, axis=0))
                    energies.append(np.mean(norms**2))
                # Requiring all three blocks prevents score-driven subgroup deletion.
                if any(v is None for v in means):
                    continue
                for i in range(3):
                    for k in range(i):
                        rows.append(dict(branch=c, angle_bin=b, channel=j, block1=k, block2=i,
                                         phase=float(means[i]@means[k]),
                                         energy=float(raw_means[i]@raw_means[k]/((energies[i]+energies[k])/2))))
    return rows


def frozen_svd(w, events, mask, cfg):
    result = w.copy()
    bins = (np.mod(events['psi'], 2*np.pi)/(np.pi/3)).astype(int)
    cal = interval(events['t'], cfg['split_seconds']['calibration'])
    fitted = 0
    for c in [0, 1]:
        for b in range(6):
            train = cal & mask & (events['branch'] == c) & (bins == b)
            use = (events['branch'] == c) & (bins == b)
            if train.sum() < cfg['svd']['minimum_events_per_fit_group']:
                continue
            for j in range(3):
                xx = w[train, j].astype(float)
                xx -= xx.mean(axis=1, keepdims=True)
                _, _, v = np.linalg.svd(xx, full_matrices=False)
                v = v[:cfg['svd']['initial_rank']]
                yy = w[use, j].astype(float)
                yy -= yy.mean(axis=1, keepdims=True)
                result[use, j] = (yy@v.T)@v
                fitted += 1
    return result, fitted


def family_spectrum(y, fs, span, motion, zr=84, zp=31):
    t = np.arange(len(y))/fs+span[0]
    # Envelope fusion is kept distinct from signed vibration fusion.
    analytic = signal.hilbert(y, axis=0)
    envelope = np.abs(analytic).mean(axis=1)
    psi = np.interp(t, motion['time'], motion['phase'])/zr
    revolutions = (psi-psi[0])/(2*np.pi)
    count = max(256, int((revolutions[-1]-revolutions[0])*4096))
    grid = np.linspace(revolutions[0], revolutions[-1], count, endpoint=False)
    e = np.interp(grid, revolutions, envelope)
    e -= e.mean()
    window = np.hanning(len(e))
    amplitude = 2*np.abs(np.fft.rfft(e*window))/window.sum()
    order = np.fft.rfftfreq(len(e), (grid[1]-grid[0]))
    scores = []
    for factor in [1, 2]:
        target = zr/zp*factor*np.arange(1, 5)
        line = np.array([amplitude[np.abs(order-f) <= .18].max(initial=0) for f in target])
        use = (order >= 1) & (order <= 30)
        for f in target:
            use &= np.abs(order-f) > .45
        noise = np.median(amplitude[use])
        scores.append(float(20*np.log10(max(np.sqrt(np.mean(line**2)), 1e-20)/max(noise, 1e-20))))
    use = order <= 30
    return order[use], amplitude[use], scores
