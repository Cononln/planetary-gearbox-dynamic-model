"""Sun-fault trial of signed-waveform channel alignment and full-tidal TSA."""
from pathlib import Path
from datetime import datetime
import argparse
import os
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0, str(ROOT/'src_py'))
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '2')
import numpy as np
from scipy import io, signal
from event_alignment import rational_resample, estimate_motion, sample_sinc, pulse_edges
from pearson_scope import pearson_lag, matrix_tsa
from signed_waveform import signed_summary, peak_preserving_indices
from tidal_waveform import recurrence, angular_rows, split_row_correlation
from test_signed_waveform_route import dump_json, dump_csv, read_csv, load_npz, digest, forbidden
from plot_fault_event_alignment import plt, save, mark

METHODS = ['none', 'whole_waveform', 'tidal_waveform']
LABELS = ['Unaligned', 'Whole-span Pearson', 'Per-tidal Pearson']
COLORS = ['#858D96', '#648EA3', '#B16D52']
CN = ['未对齐', '整段波形对齐', '逐太阳轮潮汐对齐']
RANKS = [0, 1, 2]
FS = 51200
POINTS = 16384


def prepare(source, out):
    print('Reading sun record; inherited sampling rate 51200 Hz', flush=True)
    stamp = source.stat()
    loaded = io.loadmat(source, variable_names=['chanvals'])
    values = loaded['chanvals']
    assert values.ndim == 2 and values.shape[1] == 5
    size = values.shape[0]
    first, last = 30*FS, 90*FS
    assert last <= size
    high = np.array(values[first:last, :3], dtype=np.float32, order='C')
    enc = [np.array(values[first:last, j], dtype=np.float64) for j in [3, 4]]
    del values, loaded
    assert np.isfinite(high).all()
    low = rational_resample(high, FS, 5120)
    mean = low[3*5120:15*5120].mean(axis=0)
    high -= mean
    low -= mean
    cfg = dict(working_rates_hz=dict(motion=5120), split_seconds=dict(calibration=[3, 15]),
               motion=dict(mesh_search_hz=[145, 190], harmonics=[1, 2, 3, 4]))
    motion = estimate_motion(low, cfg)
    print('Vibration-derived mesh frequency', motion['fm'], flush=True)
    np.savez(out/'motion.npz', **motion)
    # Equal to the real part of the previous complex bandpass convolution.
    n = np.arange(513)-256
    taps = 2*signal.firwin(513, 3000, fs=FS, window=('kaiser', 8))*np.cos(2*np.pi*7000*n/FS)
    wave = np.empty_like(high)
    band_fraction = []
    for j in range(3):
        wave[:, j] = signal.fftconvolve(high[:, j], taps, mode='same')
        band_fraction.append(float(np.mean(wave[:, j].astype(float)**2)/np.mean(high[:, j].astype(float)**2)))
    np.save(out/'band_waveform.npy', wave)
    del high, low
    # Encoder is inspected only after the vibration motion estimate is frozen.
    encoder_stats = []
    for j, raw in enumerate(enc):
        edges = pulse_edges(raw, FS)
        if len(edges) > 10:
            pulse_rate = (len(edges)-1)/(edges[-1]-edges[0])
            rpm = pulse_rate/1024*60
            encoder_stats.append(dict(channel=j+4, edges=int(len(edges)), pulse_rate_hz=float(pulse_rate),
                                      carrier_rpm_assuming_1024_ppr=float(rpm),
                                      implied_mesh_hz=float(pulse_rate/1024*84)))
    del enc
    info = dict(source=str(source), source_samples=size, source_channels=5,
                fs_hz=FS, fs_basis='Inherited 51200 Hz from legacy MATLAB loader; no sample-rate field in MAT',
                total_duration_s=size/FS, source_window_s=[30, 90], vibration_columns_one_based=[1, 2, 3],
                encoder_posterior_only=encoder_stats, vibration_mesh_hz=float(motion['fm']),
                band_energy_fraction=band_fraction, source_bytes=stamp.st_size,
                source_mtime_ns=stamp.st_mtime_ns, source_sha256=digest(source))
    dump_json(out/'source_audit.json', info)
    return wave, motion, info


def compute(source, out):
    signal.hilbert = forbidden
    wave, motion, info = prepare(source, out)
    kin = recurrence('sun')
    edges, query = angular_rows(motion, kin['carrier_turns'], POINTS)
    count = len(edges)-1
    duration = float(np.mean(np.diff(edges)))
    assert POINTS/max(np.diff(edges)) > 22000
    print('SUN TIDAL ROWS', count, 'duration', duration, 's; points', POINTS, flush=True)
    dump_csv(out/'periods.csv', [dict(tidal=q+1, local_start_s=edges[q], local_end_s=edges[q+1],
                                     source_start_s=edges[q]+30, source_end_s=edges[q+1]+30,
                                     duration_s=edges[q+1]-edges[q]) for q in range(count)])
    shifts = {method: np.zeros((count, 3)) for method in METHODS}
    fits, curves, checks, metrics, energies = [], [], [], [], []
    for method in METHODS[1:]:
        spans = [(edges[0], edges[-1])] if method.startswith('whole') else list(zip(edges[:-1], edges[1:]))
        for q, (start, end) in enumerate(spans):
            a, b = round(start*FS), round(end*FS)
            for j in [1, 2]:
                result, lag, corr = pearson_lag(wave[a:b, 0], wave[a:b, j], 410)
                shift = result['lag_samples']/FS
                if method.startswith('whole'):
                    shifts[method][:, j] = shift
                else:
                    shifts[method][q, j] = shift
                identity = dict(method=method, tidal=0 if method.startswith('whole') else q+1, channel=j+1)
                fits.append(dict(**identity, lag_ms=float(shift*1000), pearson_zero=result['zero'],
                                 pearson_peak=result['peak'], boundary=result['boundary'], valid=result['valid']))
                # Preserve every fitted curve numerically, including ambiguous secondary peaks.
                curves.extend(dict(**identity, lag_ms=float(a/FS*1000), pearson=float(b)) for a, b in zip(lag, corr))
        print(method, 'lag fit complete', flush=True)
    saved, fixed_s1, numerical_error = {}, {}, []
    for method in METHODS:
        print(method, 'resampling', count, 'signed rows', flush=True)
        matrix = np.empty((count, 3, POINTS), dtype=np.float32)
        for j in range(3):
            samples = (query+shifts[method][:, j, None])*FS
            assert samples.min() > 16 and samples.max() < len(wave)-17
            matrix[:, j] = sample_sinc(wave[:, j], samples).astype(np.float32)
            checks.append(dict(method=method, channel=j+1, odd_even_row_tsa_correlation=split_row_correlation(matrix[:, j])))
        checks.append(dict(method=method, channel=0, odd_even_row_tsa_correlation=split_row_correlation(matrix.mean(axis=1))))
        saved[method+'_lags_ms'] = (shifts[method]*1000).astype(np.float32)
        for rank in RANKS:
            channels = []
            for j in range(3):
                mean, singular = matrix_tsa(matrix[:, j], rank)
                channels.append(mean)
                if rank == 1:
                    energies.append(dict(method=method, channel=j+1, rank1_fraction=float(singular[0]**2/np.sum(singular**2))))
            tsa = np.asarray(channels, dtype=np.float32)
            if method == 'none':
                fixed_s1[rank] = tsa[0].copy()
            else:
                assert np.array_equal(tsa[0], fixed_s1[rank])
            fused, frequency, amplitude, row = signed_summary(tsa, duration, fault_bins=kin['fault_cycles'],
                                                              carrier_bins=kin['carrier_turns'], mesh_bins=kin['mesh_cycles'])
            key = f'{method}_rank{rank}'
            saved[key+'_channels'] = tsa
            saved[key+'_fused'] = fused.astype(np.float32)
            saved[key+'_amplitude'] = amplitude.astype(np.float32)
            assert np.isfinite(fused).all() and 0 <= row['channel_coherent_energy_fraction'] <= 1+1e-12
            parseval = .5*np.sum(amplitude[1:-1]**2)+amplitude[0]**2+amplitude[-1]**2
            numerical_error.append(abs(parseval-np.mean(fused*fused))/max(np.mean(fused*fused), 1e-30))
            metrics.append(dict(method=method, rank=rank, tidal_rows=count, **row))
        del matrix
        print(method, 'SVD/TSA done', flush=True)
    saved['frequency_hz'] = frequency
    saved['equivalent_time_s'] = np.arange(POINTS)/POINTS*duration
    saved['duration_s'] = np.asarray(duration)
    np.savez_compressed(out/'Sun050_signed_results.npz', **saved)
    dump_csv(out/'metrics.csv', metrics)
    dump_csv(out/'lag_fits.csv', fits)
    dump_csv(out/'lag_curves.csv', curves)
    dump_csv(out/'row_repeatability.csv', checks)
    dump_csv(out/'svd_energy.csv', energies)
    assert max(numerical_error) < 1e-10
    assert digest(source) == info['source_sha256']
    dump_json(out/'run_manifest.json', dict(
        complete=True, source=info, fault='sun', kinematics=kin, geometry_status='inherited 21/31/84, audit pending',
        rows=count, row_duration_s=duration, points_per_row=POINTS,
        analyzed_source_span_s=[float(edges[0]+30), float(edges[-1]+30)],
        matching='Signed waveform Pearson, S2/S3 to fixed S1, one lag per chosen span',
        processing='4-10 kHz signed vibration -> angular rows -> per-channel SVD/TSA -> signed channel mean -> direct FFT',
        psi='not used in this pilot; original PSI requires modal envelope',
        original_waveform_envelope=False, individual_event_registration=False, extra_row_registration=False,
        motion='Vibration multi-harmonic phase tracking; complex magnitudes weight phase increments only, not TSA data',
        ranks=RANKS, evidence='Descriptive single historical record; row count is not independent record count',
        code_sha256={str(p.relative_to(ROOT)): digest(p) for p in
                     [Path(__file__), ROOT/'src_py'/'tidal_waveform.py', ROOT/'src_py'/'signed_waveform.py', ROOT/'src_py'/'pearson_scope.py']},
    ))
    dump_json(out/'numerical_audit.json', dict(passed=True, full_tidal_rows=count, s1_unchanged=True,
                                              maximum_relative_parseval_error=max(numerical_error),
                                              all_lag_fits_valid=all(r['valid'] for r in fits),
                                              boundary_hits=sum(r['boundary'] for r in fits),
                                              source_unchanged=True, no_analytic_envelope_call=True))


def plot(out):
    data = load_npz(out/'Sun050_signed_results.npz')
    t, f = data['equivalent_time_s'], data['frequency_hz']
    duration = float(data['duration_s'])
    rows = read_csv(out/'periods.csv')
    n = len(rows)
    baseline = data['none_rank0_channels'][0]
    interior = np.flatnonzero((t >= .005) & (t <= duration-.005))
    center = float(t[interior[np.argmax(abs(baseline[interior]))]])
    zoom = np.flatnonzero(abs(t-center) <= .002)
    qa, source = [], []
    fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 126/25.4))
    fig.subplots_adjust(left=.11, right=.98, top=.87, bottom=.15, hspace=.57, wspace=.36)
    for r, rank in enumerate([0, 1]):
        for c, ax in enumerate(axes[r]):
            for method, label, color in zip(METHODS, LABELS, COLORS):
                y = data[f'{method}_rank{rank}_fused']
                ix = peak_preserving_indices(y, 0, len(y), 1000) if c == 0 else zoom
                x = t[ix] if c == 0 else (t[ix]-center)*1000
                ax.plot(x, y[ix], color=color, lw=.65 if c == 0 else .9, label=label)
                source.extend(dict(panel=f'wave_{r}_{c}', method=method, rank=rank, x=float(a), value=float(b)) for a, b in zip(x, y[ix]))
            ax.axhline(0, color='.65', lw=.4)
            ax.set(xlabel='Equivalent time (s)' if c == 0 else 'Equivalent local time (ms)',
                   ylabel='Signed acquisition amplitude', xlim=(0, duration) if c == 0 else (-2, 2))
            mark(ax, chr(97+2*r+c), ('Direct TSA' if rank == 0 else 'Rank-1 SVD + TSA')+(' | full tidal row' if c == 0 else ' | detail'))
    limit = max(abs(x) for ax in axes.flat for x in ax.get_ylim())
    for ax in axes.flat:
        ax.set_ylim(-limit, limit)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.55, .985), ncol=3)
    fig.text(.11, .035, f'Sun050, 600 rpm, L00: {n} tidal rows, {duration:.4f} s each; signed three-channel fusion.\nShared amplitude scale. Zoom fixed from unshifted S1; not confirmed contact timing. No artificial repetition.', fontsize=6, color='.3')
    save(fig, out, '01_sun_signed_tsa', qa)

    ba = data['none_rank0_amplitude']
    band = np.flatnonzero((f >= 4000) & (f <= 10000))
    spectral_center = float(f[band[np.argmax(ba[band])]])
    fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 125/25.4))
    fig.subplots_adjust(left=.11, right=.98, top=.87, bottom=.15, hspace=.57, wspace=.36)
    for r, rank in enumerate([0, 1]):
        for c, ax in enumerate(axes[r]):
            for method, label, color in zip(METHODS, LABELS, COLORS):
                y = data[f'{method}_rank{rank}_amplitude']
                ix = np.flatnonzero(f <= 12000) if c == 0 else np.flatnonzero(abs(f-spectral_center) <= 150)
                x = f[ix]/1000 if c == 0 else f[ix]
                ax.plot(x, y[ix], color=color, lw=.7, label=label)
                source.extend(dict(panel=f'spectrum_{r}_{c}', method=method, rank=rank, x=float(a), value=float(b)) for a, b in zip(x, y[ix]))
            ax.set(xlabel='Equivalent frequency (kHz)' if c == 0 else 'Equivalent frequency (Hz)',
                   ylabel='Direct Fourier amplitude', xlim=(0, 12) if c == 0 else (spectral_center-150, spectral_center+150))
            mark(ax, chr(97+2*r+c), ('Direct TSA' if rank == 0 else 'Rank-1 SVD + TSA')+(' | spectrum' if c == 0 else ' | detail'))
    upper = max(ax.get_ylim()[1] for ax in axes.flat)
    for ax in axes.flat:
        ax.set_ylim(0, upper)
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.55, .985), ncol=3)
    fig.text(.11, .035, f'Direct spectrum of signed fused TSA. Fourier grid {1/duration:.3f} Hz; sun-fault spacing {12/duration:.3f} Hz.\nMore rows improve averaging, not Fourier resolution. Band fixed at 4-10 kHz; all y axes start at zero.', fontsize=6, color='.3')
    save(fig, out, '02_sun_direct_spectra', qa)

    fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 87/25.4))
    fig.subplots_adjust(left=.105, right=.98, top=.82, bottom=.24, wspace=.35)
    for ax, channel in zip(axes, [2, 3]):
        for method, label, color in zip(METHODS[1:], LABELS[1:], COLORS[1:]):
            y = data[method+'_lags_ms'][:, channel-1]
            ax.plot(np.arange(1, n+1), y, color=color, lw=.75, label=label)
        ax.axhline(0, color='.65', lw=.4)
        ax.set(xlabel='Sun tidal row', ylabel='Waveform lag (ms)', xlim=(1, n))
        mark(ax, chr(95+channel), f'S{channel} relative to fixed S1')
    hh, ll = axes[0].get_legend_handles_labels()
    fig.legend(hh, ll, loc='upper center', bbox_to_anchor=(.55, .98), ncol=2)
    fig.text(.105, .035, 'One scalar lag per full tidal row, not individual-pulse alignment. S1 remains fixed.\nLarge jumps may indicate alternative correlation peaks; fitted similarity is not independent phase truth.', fontsize=6, color='.3')
    save(fig, out, '03_sun_tidal_lags', qa)
    dump_csv(out/'figure_source_data.csv', source)
    dump_json(out/'figure_qa.json', dict(backend='Python', figures=qa, visual_inspection='pending',
                                        waveform_zoom_s=center, spectral_zoom_hz=spectral_center))


def report(out):
    import json
    m = json.loads((out/'run_manifest.json').read_text(encoding='utf-8'))
    rows = read_csv(out/'metrics.csv')
    checks = read_csv(out/'row_repeatability.csv')
    text = ['# 太阳轮故障有符号波形试验', '',
            f'使用太阳轮 50% 故障、600 rpm、L00 三通道记录；本轮未处理行星轮故障信号。',
            f'按既有 51,200 Hz 解释，原记录长 {m["source"]["total_duration_s"]:.3f} s；MAT 内没有采样率字段，因此该值属于原脚本沿用假设。提取原记录 30-90 s，最终分析范围为 {m["analyzed_source_span_s"][0]:.5f}-{m["analyzed_source_span_s"][1]:.5f} s。', '',
            '## TSA 周期与实际平均', '',
            f'- 按太阳轮21、齿圈84、三行星的当前配置：太阳轮潮汐周期为1个行星架转动周期，包含12次太阳轮故障啮合事件和84个啮合周期。该定义不要求每个健康行星齿的身份同时复现。',
            f'- 每行约 {m["row_duration_s"]:.6f} s，重采样为 {m["points_per_row"]} 点；共平均 {m["rows"]} 行。每通道矩阵为 {m["rows"]} × {m["points_per_row"]}。',
            '- 对齐仍是每行内 S2/S3 对 S1，整段法仅各一个时移；没有额外将不同潮汐行或逐故障事件相互配准。',
            '- 保留带正负号的4-10 kHz振动，依次执行波形相关匹配、角域重采样、每通道SVD/直接TSA、三通道有符号融合、直接FFT。原稿PSI依赖包络，本轮仍不使用，固定秩0/1/2为预试验对照。',
            f'- 输出仍只有一个约 {m["row_duration_s"]:.4f} s 的周期，不能将其人为重复当作1 s实测波形。频谱网格约 {1/m["row_duration_s"]:.3f} Hz；平均次数增加不等于频率分辨率提高。', '',
            '## 结果', '', '| SVD秩 | 方法 | 融合RMS | 通道相干能量比例 | 故障间隔梳状谱/背景 dB |', '|---:|---|---:|---:|---:|']
    for rank in RANKS:
        for method, name in zip(METHODS, CN):
            r = next(v for v in rows if v['method'] == method and int(v['rank']) == rank)
            text.append(f'| {rank} | {name} | {float(r["signed_fused_rms"]):.6f} | {float(r["channel_coherent_energy_fraction"]):.4f} | {float(r["fault_spaced_comb_background_db"]):.3f} |')
    text += ['', '## 奇偶潮汐行平均波形的一致性', '', '| 方法 | S1 | S2 | S3 | 融合 |', '|---|---:|---:|---:|---:|']
    for method, name in zip(METHODS, CN):
        val = [float(next(r['odd_even_row_tsa_correlation'] for r in checks if r['method'] == method and int(r['channel']) == j)) for j in [1, 2, 3, 0]]
        text.append('| '+name+' | '+' | '.join(f'{v:.4f}' for v in val)+' |')
    text += ['', '## 解释边界', '',
             '- 相干能量比例衡量通道是否更容易相加，不是故障选择性证据。奇偶行一致性为单条记录的描述性重复性检查，匹配仍使用对应数据，不是独立盲测。',
             '- 故障间隔梳状谱指标在4-10 kHz中统计每12个Fourier bin的目标及左右各1 bin，排除啮合线附近；匹配背景为目标左右5-10 bin。它不是包络谱，也不是物理SNR。太阳轮与上一轮行星轮的目标周期、网格和数据不同，数值不能直接横向判优。',
             '- 正常啮合也可能同步保留。幅值下降既可能是干扰衰减，也可能是目标被抵消；本轮没有真实故障接触时刻标签。',
             '- 频带沿用4-10 kHz，未根据太阳轮结果重新挑选；三通道极性和幅值保持原样。当前齿数几何一致性及灵敏度标定仍待核实。',
             '- 来源元数据与编码器事后核查见 source_audit.json；完整波形、谱、时移和全部相关曲线均已保留。', '']
    (out/'太阳轮无包络试验报告.md').write_text('\n'.join(text), encoding='utf-8')
    for r in rows:
        print('RESULT', r['method'], 'rank', r['rank'], 'coherence', r['channel_coherent_energy_fraction'],
              'comb_dB', r['fault_spaced_comb_background_db'], flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True,
                        help='Local five-column MAT record (not stored in Git)')
    parser.add_argument('--plot-only', type=Path)
    args = parser.parse_args()
    if args.plot_only:
        out = args.plot_only.resolve()
    else:
        out = ROOT/'results'/('sun_signed_waveform_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
        out.mkdir()
        print('OUTPUT', out, flush=True)
        compute(args.source.resolve(), out)
    report(out)
    plot(out)
    print('COMPLETE', out, flush=True)


if __name__ == '__main__':
    main()
