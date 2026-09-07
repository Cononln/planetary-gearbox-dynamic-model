"""Envelope-event timing pilot: measured candidate selection with frozen priors."""
from __future__ import annotations
import numpy as np
from scipy import signal, ndimage
from event_alignment import interval, path_design, extract_events


def envelope_response(wave, fs):
    out = np.empty(wave.shape, np.float32)
    for j in range(3):
        out[:, j] = ndimage.gaussian_filter1d(np.abs(signal.hilbert(wave[:, j])), .00012*fs)
    return out


def onset_markers(env, fs, events):
    """A fixed strongest prominent-packet marker; never use a predicted lag here."""
    onset = np.full((len(events['t']), 3), np.nan)
    quality = np.zeros_like(onset)
    for n, t0 in enumerate(events['t']):
        a, b = round((t0-.012)*fs), round((t0+.012)*fs)
        for j in range(3):
            y = env[a:b, j].astype(float)
            baseline = np.percentile(y, 25)
            low = y[y <= np.percentile(y, 60)]
            scale = max(1.4826*np.median(np.abs(low-np.median(low))), .05*baseline, 1e-12)
            peaks, props = signal.find_peaks(y, prominence=3*scale, distance=max(1, round(.0008*fs)))
            valid = (peaks+a >= (t0-.008)*fs) & (peaks+a <= (t0+.008)*fs)
            if not valid.any():
                continue
            peaks, prominence = peaks[valid], props['prominences'][valid]
            k = int(np.argmax(prominence)); p = peaks[k]
            level = baseline+.2*(y[p]-baseline)
            left = max(0, p-round(.003*fs))
            crossings = np.flatnonzero((y[left:p] <= level) & (y[left+1:p+1] > level))+left
            if not len(crossings):
                continue
            q = crossings[-1]
            tx = (a+q+(level-y[q])/max(y[q+1]-y[q], 1e-12))/fs-t0
            if abs(tx) <= .008:
                onset[n, j] = tx
                quality[n, j] = prominence[k]/scale
    return dict(offset=onset, quality=quality)


def arrival_templates(env, fs, events, markers, calibration):
    offsets = np.nan_to_num(markers['offset'])
    windows, _ = extract_events(env, fs, events, offsets)
    bins = (np.mod(events['psi'], 2*np.pi)/(np.pi/3)).astype(int)
    templates, audit = {}, []
    for c in [0, 1]:
        for b in range(6):
            for j in range(3):
                base = interval(events['t'], calibration) & (events['branch'] == c) & np.isfinite(markers['offset'][:, j])
                use = base & (bins == b)
                fallback = use.sum() < 4
                if fallback:
                    use = base
                if use.sum() < 4:
                    continue
                yy = windows[use, j].astype(float)
                yy -= np.percentile(yy, 20, axis=1, keepdims=True)
                yy /= np.maximum(np.linalg.norm(yy, axis=1, keepdims=True), 1e-12)
                template = np.median(yy, axis=0)
                template -= template.mean()
                templates[(c, b, j)] = template
                audit.append(dict(branch=c, angle_bin=b, channel=j, calibration_events=int(use.sum()), branch_fallback=fallback))
    return templates, audit


def envelope_candidates(env, fs, events, templates, support=None):
    n = len(events['t']); kmax = 4
    offsets = np.full((n, 3, kmax), np.nan)
    scores = np.full_like(offsets, -np.inf)
    width, search = round(.020*fs), round(.008*fs)
    bins = (np.mod(events['psi'], 2*np.pi)/(np.pi/3)).astype(int)
    for i, (t0, c, b) in enumerate(zip(events['t'], events['branch'], bins)):
        start = round((t0-.004)*fs)
        for j in range(3):
            if support is not None and not support[i, j]:
                continue
            template = templates.get((int(c), int(b), j))
            if template is None or np.linalg.norm(template) < 1e-15:
                continue
            y = env[start-search:start+width+search, j].astype(float)
            ss = np.r_[0, np.cumsum(y)]; sq = np.r_[0, np.cumsum(y*y)]
            variance = sq[width:]-sq[:-width]-(ss[width:]-ss[:-width])**2/width
            corr = signal.correlate(y, template, mode='valid', method='fft')
            corr /= np.maximum(np.sqrt(np.maximum(variance, 0))*np.linalg.norm(template), 1e-15)
            peaks, _ = signal.find_peaks(np.r_[-np.inf, corr, -np.inf], height=.25,
                                         distance=max(1, round(.0008*fs)))
            peaks -= 1
            peaks = peaks[np.argsort(corr[peaks])[::-1]][:kmax]
            for k, p in enumerate(peaks):
                sub = 0.
                if 0 < p < len(corr)-1:
                    den = corr[p-1]-2*corr[p]+corr[p+1]
                    sub = np.clip(.5*(corr[p-1]-corr[p+1])/min(den, -1e-12), -.5, .5)
                offsets[i, j, k] = np.clip((p-search+sub)/fs+start/fs-(t0-.004), -.008, .008)
                scores[i, j, k] = corr[p]
    return dict(offsets=offsets, scores=scores)


def local_result(candidates):
    valid = np.isfinite(candidates['offsets'][:, :, 0])
    return dict(shift=np.nan_to_num(candidates['offsets'][:, :, 0]), channel_valid=valid,
                accepted=valid.sum(axis=1) >= 2)


def fit_relative_path(events, candidates, calibration):
    use = interval(events['t'], calibration)
    basis = path_design(events['psi'][use], events['branch'][use])
    y = candidates['offsets'][use, :, 0]*1000
    q = np.clip(candidates['scores'][use, :, 0], 0, 1)**2
    valid = np.isfinite(y) & (np.isfinite(y).sum(axis=1)[:, None] >= 2)
    weight = np.where(valid, q, 0)
    yy = np.nan_to_num(y)
    beta = np.zeros(20)
    regularizer = np.tile([.001, .001, .03, .03, .12, .12, .15, .15, .15, .15], 2)
    for _ in range(30):
        d = np.sum(weight*(yy-basis@beta), axis=1)/np.maximum(weight.sum(axis=1), 1e-12)
        residual = yy-d[:, None]-basis@beta
        robust = np.minimum(1., .7/np.maximum(np.abs(residual), 1e-12))*weight
        x = basis.reshape(-1, 20); w = robust.ravel()
        beta = np.linalg.solve(x.T@(w[:, None]*x)+np.diag(regularizer),
                               x.T@(w*(yy-d[:, None]).ravel()))
    grid = path_design(np.tile(np.linspace(0, 2*np.pi, 145), 2), np.repeat([0, 1], 145))
    maximum = np.max(np.abs(grid@beta))
    if maximum > 4:
        beta *= 4/maximum
    # Simpler calibration-only channel-offset baseline for held-out prediction.
    const = np.zeros((2, 3))
    for c in [0, 1]:
        rows = (events['branch'][use] == c) & (valid.sum(axis=1) == 3)
        if rows.any():
            relative = yy[rows]-yy[rows].mean(axis=1, keepdims=True)
            const[c] = np.median(relative, axis=0)
            const[c] -= const[c].mean()
    return dict(coeff=beta, constant_ms=const, calibration_supported_events=int((valid.sum(axis=1) >= 2).sum()))


def huber(x):
    a = np.abs(x)
    return np.where(a < 1, .5*x*x, a-.5)


def joint_result(events, cand, model, spans, strength=.5, excluded=None):
    path = path_design(events['psi'], events['branch'])@np.asarray(model['coeff'])
    n = len(events['t']); states = np.linspace(-8, 8, 81)
    shift = np.zeros((n, 3)); channel_valid = np.zeros((n, 3), bool)
    common = np.zeros(n); assignment = np.full((n, 3), -1, int)
    prediction = np.zeros((n, 3))
    for span in spans:
        ids = np.flatnonzero(interval(events['t'], span))
        if not len(ids):
            continue
        lag = cand['offsets'][ids]*1000
        score = cand['scores'][ids].copy()
        if excluded is not None:
            lag = lag.copy(); lag[:, excluded] = np.nan; score[:, excluded] = -np.inf
        best = np.max(score, axis=2, keepdims=True)
        # Invalid-candidate arithmetic is explicitly masked below.
        finite = np.isfinite(lag) & np.isfinite(score)
        likelihood = (np.where(np.isfinite(best), best, 0)-np.where(finite, score, 0))/.12
        likelihood += -.5*np.log(np.maximum(np.where(finite, score, 0), 1e-6))
        residual = (lag[:, None]-states[None, :, None, None]-path[ids, None, :, None])/.7
        cost = likelihood[:, None]+strength*huber(np.nan_to_num(residual))
        cost = np.where(finite[:, None], cost, np.inf)
        chosen = np.argmin(cost, axis=3)
        minimum = np.min(cost, axis=3)
        valid = minimum < 2.5
        unary = np.minimum(minimum, 2.5).sum(axis=2)
        unary += 5*(valid.sum(axis=2) < 2)
        legal = np.all(np.abs(states[None, :, None]+path[ids, None]) <= 8.00001, axis=2)
        unary = np.where(legal, unary, 1e6)
        back = np.zeros((len(ids), len(states)), np.int32)
        value = unary[0].copy()
        ref = np.median(np.diff(events['t'][ids])) if len(ids) > 1 else .1
        for i in range(1, len(ids)):
            scale = max(.5, np.sqrt((events['t'][ids[i]]-events['t'][ids[i-1]])/ref))
            transition = .08*np.minimum(huber((states[:, None]-states[None, :])/scale), 8)
            total = value[:, None]+transition
            back[i] = np.argmin(total, axis=0)
            value = unary[i]+np.min(total, axis=0)
        state = np.zeros(len(ids), int); state[-1] = np.argmin(value)
        for i in range(len(ids)-1, 0, -1):
            state[i-1] = back[i, state[i]]
        for i, index in enumerate(ids):
            g = state[i]; common[index] = states[g]/1000
            prediction[index] = (states[g]+path[index])/1000
            for j in range(3):
                if valid[i, g, j]:
                    k = chosen[i, g, j]
                    shift[index, j] = lag[i, j, k]/1000
                    channel_valid[index, j] = True; assignment[index, j] = k
    accepted = channel_valid.sum(axis=1) >= 2
    shift[~accepted] = 0
    return dict(shift=shift, accepted=accepted, channel_valid=channel_valid, common=common,
                predicted=prediction, path=path/1000, assignment=assignment)


def other_channel_baseline(events, cand, model, excluded):
    offsets = cand['offsets'][:, :, 0]*1000
    const = np.asarray(model['constant_ms'])[events['branch']]
    others = [j for j in range(3) if j != excluded]
    valid = np.all(np.isfinite(offsets[:, others]), axis=1)
    common = np.mean(np.nan_to_num(offsets[:, others]-const[:, others]), axis=1)
    prediction = np.clip((common+const[:, excluded])/1000, -.008, .008)
    return prediction, valid


def cross_channel_marker_error(markers, shifts, event_mask):
    residual = markers['offset']-shifts
    values = []
    for n in np.flatnonzero(event_mask):
        for j in range(3):
            for k in range(j):
                if np.isfinite(residual[n, j]) and np.isfinite(residual[n, k]):
                    values.append(abs(residual[n, j]-residual[n, k])*1000)
    return float(np.mean(values)) if values else np.nan, len(values)
