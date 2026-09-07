"""Recompute a signed-waveform scope pilot and export auditable comparisons."""
from pathlib import Path
from datetime import datetime
import argparse
import csv
import hashlib
import json
import os
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0, str(ROOT/'src_py'))
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '2')
import numpy as np
from scipy import signal
from event_alignment import sample_sinc
from pearson_scope import pearson_lag, matrix_tsa
from signed_waveform import signed_summary, peak_preserving_indices
from plot_fault_event_alignment import plt, save, mark

METHODS = ['none', 'whole_waveform', 'tidal_waveform']
LABELS = ['Unaligned', 'Whole-span Pearson', 'Per-tidal Pearson']
COLORS = ['#858D96', '#648EA3', '#B16D52']
CN = ['未对齐', '整段波形对齐', '逐潮汐波形对齐']
POINTS = 524288
RANKS = [0, 1, 2]


def dump_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')


def dump_csv(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def load_npz(path):
    with np.load(path) as data:
        return {k: data[k] for k in data.files}


def digest(path):
    sha = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            sha.update(block)
    return sha.hexdigest()


def forbidden(*args, **kwargs):
    raise AssertionError('Analytic-envelope calculation is disabled in this pilot')


def find_periods(motion, ring, planet):
    carrier = motion['phase']/(2*np.pi*ring)
    repeat = planet//int(np.gcd(planet, ring))
    origin = np.interp(3., motion['time'], carrier)
    count = int((np.interp(59., motion['time'], carrier)-origin)//repeat)
    assert count == 3
    phase_edges = origin+np.arange(count+1)*repeat
    edges = np.interp(phase_edges, carrier, motion['time'])
    within = np.arange(POINTS)/POINTS*repeat
    query = np.interp((phase_edges[:-1, None]+within).ravel(), carrier, motion['time'])
    return edges, query.reshape(count, POINTS), repeat


def compute(source, prior, out):
    signal.hilbert = forbidden
    cfg = json.loads((source/'runtime_config.json').read_text(encoding='utf-8'))
    fs = cfg['working_rates_hz']['impulse']
    ring, planet = cfg['gear_teeth']['ring'], cfg['gear_teeth']['planet']
    assert (ring, planet, fs) == (84, 31, 51200)
    half = round(.008*fs)
    metrics, periods, lags, curves, energy = [], [], [], [], []
    hashes, cache_errors, summary_errors = {}, [], []
    for condition in ['BL', 'PF50']:
        wave_path = source/condition/'band_2'/'waveform.npy'
        motion_path = source/condition/'motion.npz'
        hashes[str(wave_path)] = digest(wave_path)
        hashes[str(motion_path)] = digest(motion_path)
        wave = np.load(wave_path, mmap_mode='r')
        motion = load_npz(motion_path)
        edges, query, repeat = find_periods(motion, ring, planet)
        duration = float(np.mean(np.diff(edges)))
        assert POINTS/np.max(np.diff(edges)) > 22000
        for q in range(3):
            periods.append(dict(condition=condition, tidal=q+1, start_s=edges[q],
                                end_s=edges[q+1], duration_s=edges[q+1]-edges[q]))
        shifts = {method: np.zeros((3, 3)) for method in METHODS}
        for method in METHODS[1:]:
            spans = [(edges[0], edges[-1])] if method.startswith('whole') else list(zip(edges[:-1], edges[1:]))
            for q, (start, end) in enumerate(spans):
                a, b = round(start*fs), round(end*fs)
                for j in [1, 2]:
                    result, grid, corr = pearson_lag(wave[a:b, 0], wave[a:b, j], half)
                    value = result['lag_samples']/fs
                    if method.startswith('whole'):
                        shifts[method][:, j] = value
                    else:
                        shifts[method][q, j] = value
                    identity = dict(condition=condition, method=method,
                                    tidal=0 if method.startswith('whole') else q+1, channel=j+1)
                    lags.append(dict(**identity, lag_ms=1000*value, pearson_zero=result['zero'],
                                     pearson_peak=result['peak'], boundary=result['boundary'], valid=result['valid']))
                    curves.extend(dict(**identity, lag_ms=float(g/fs*1000), pearson=float(c))
                                  for g, c in zip(grid, corr))
            print(condition, method, 'signed lags (ms)', (shifts[method]*1000).tolist(), flush=True)
        data = {}
        fixed_s1 = {}
        for method in METHODS:
            print(condition, method, 'resampling signed rows', flush=True)
            matrix = np.empty((3, 3, POINTS), np.float32)
            for j in range(3):
                samples = (query+shifts[method][:, j, None])*fs
                assert samples.min() > 16 and samples.max() < len(wave)-17
                matrix[:, j] = sample_sinc(wave[:, j], samples).astype(np.float32)
            cached_path = prior/condition/f'{method}_tsa.npz'
            hashes[str(cached_path)] = digest(cached_path)
            with np.load(cached_path) as cache:
                for rank in RANKS:
                    channels = []
                    for j in range(3):
                        mean, singular = matrix_tsa(matrix[:, j], rank)
                        channels.append(mean)
                        if rank == 1:
                            energy.append(dict(condition=condition, method=method, channel=j+1,
                                               rank1_fraction=float(singular[0]**2/np.sum(singular**2))))
                    tsa = np.asarray(channels, np.float32)
                    cached = cache[f'waveform_rank{rank}']
                    error = float(np.max(abs(tsa-cached)))
                    cache_errors.append(error)
                    assert error < 1e-6*max(1., float(np.max(abs(cached))))
                    if method == 'none':
                        fixed_s1[rank] = tsa[0].copy()
                    else:
                        assert np.array_equal(tsa[0], fixed_s1[rank])
                    fused, frequency, amplitude, row = signed_summary(tsa, duration)
                    assert np.isfinite(fused).all() and np.isfinite(amplitude).all()
                    assert 0 <= row['channel_coherent_energy_fraction'] <= 1+1e-12
                    # Parseval for one-sided Fourier amplitudes, excluding DC/Nyquist.
                    spectral = .5*np.sum(amplitude[1:-1]**2)+amplitude[0]**2+amplitude[-1]**2
                    parseval = abs(spectral-np.mean(fused*fused))/max(np.mean(fused*fused), 1e-30)
                    summary_errors.append(parseval)
                    assert parseval < 1e-10
                    data[f'{method}_rank{rank}_channels'] = tsa
                    data[f'{method}_rank{rank}_fused'] = fused.astype(np.float32)
                    data[f'{method}_rank{rank}_amplitude'] = amplitude.astype(np.float32)
                    metrics.append(dict(condition=condition, method=method, rank=rank, **row))
            del matrix
            print(condition, method, 'signed SVD/TSA and direct FFT complete', flush=True)
        data['frequency_hz'] = frequency
        data['equivalent_time_s'] = np.arange(POINTS)/POINTS*duration
        data['duration_s'] = np.asarray(duration)
        np.savez_compressed(out/f'{condition}_signed_results.npz', **data)
        del data, wave
    for path, expected in hashes.items():
        assert digest(Path(path)) == expected
    dump_csv(out/'metrics.csv', metrics)
    dump_csv(out/'periods.csv', periods)
    dump_csv(out/'lag_fits.csv', lags)
    dump_csv(out/'lag_curves.csv', curves)
    dump_csv(out/'svd_energy.csv', energy)
    dump_json(out/'run_manifest.json', dict(
        complete=True, source=str(source), historical_signed_cache=str(prior),
        alignment='Signed waveform Pearson, fixed S1, S2/S3 to S1; whole span or one lag per full tidal row',
        processing='Signed angular rows -> per-channel optional rank-1/2 SVD -> row TSA -> signed channel mean -> direct FFT',
        envelope_calls=0, runtime_analytic_envelope_guard=True, psi='not used; original PSI depends on envelope',
        source_fs_hz=fs, band_hz=[4000, 10000], points_per_row=POINTS, rows=3,
        period='31 carrier turns; inherited planet31/ring84, geometry provisional',
        inter_row_extra_alignment=False, individual_fault_event_alignment=False,
        frequency='Equivalent angular Fourier frequency; not low-frequency envelope spectrum',
        evidence='One historical record per state, three dependent rows, in-sample descriptive pilot',
        source_hashes_unchanged=hashes,
        code_hashes={str(p.relative_to(ROOT)): digest(p) for p in
                     [Path(__file__), ROOT/'src_py'/'signed_waveform.py', ROOT/'src_py'/'pearson_scope.py', ROOT/'tests'/'test_signed_waveform.py']},
    ))
    dump_json(out/'numerical_audit.json', dict(
        passed=True, matching_and_output_are_signed=True, s1_unchanged=True,
        max_difference_from_historical_signed_tsa=max(cache_errors),
        max_relative_parseval_error=max(summary_errors),
        lag_boundary_hits=sum(row['boundary'] for row in lags),
        all_fits_valid=all(row['valid'] for row in lags), inputs_unchanged=True,
        unit_tests='signed cancellation; amplitude normalization; lag sign; SVD identity; display extrema; no analytic-envelope call',
    ))
    return metrics


def plot(out):
    data = load_npz(out/'PF50_signed_results.npz')
    time, freq = data['equivalent_time_s'], data['frequency_hz']
    baseline = data['none_rank0_channels'][0]
    one_second = np.flatnonzero(time < 1.)
    center_index = one_second[np.argmax(abs(baseline[one_second]))]
    center = float(time[center_index])
    zoom_time = np.flatnonzero((time >= center-.002) & (time <= center+.002))
    display, qa = [], []
    fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 127/25.4))
    fig.subplots_adjust(left=.11, right=.98, top=.87, bottom=.15, hspace=.56, wspace=.38)
    for r, rank in enumerate([0, 1]):
        for c, ax in enumerate(axes[r]):
            for method, label, color in zip(METHODS, LABELS, COLORS):
                y = data[f'{method}_rank{rank}_fused']
                ix = peak_preserving_indices(y, 0, len(one_second), 1200) if c == 0 else zoom_time
                x = time[ix] if c == 0 else (time[ix]-center)*1000
                ax.plot(x, y[ix], color=color, lw=.55 if c == 0 else .85, alpha=.85, label=label)
                display.extend(dict(panel=f'{r},{c}', method=method, rank=rank, x=float(a), amplitude=float(b)) for a, b in zip(x, y[ix]))
            ax.axhline(0, color='.65', lw=.4, zorder=0)
            ax.set(xlabel='Equivalent time (s)' if c == 0 else 'Equivalent local time (ms)',
                   ylabel='Signed acquisition amplitude', xlim=(0, 1) if c == 0 else (-2, 2))
            mark(ax, chr(97+2*r+c), ('Direct TSA' if rank == 0 else 'Rank-1 SVD + TSA')+(' | 1 s' if c == 0 else ' | 4 ms'))
    limit = max(abs(a) for ax in axes.flat for a in ax.get_ylim())
    for ax in axes.flat:
        ax.set_ylim(-limit, limit)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.55, .985), ncol=3)
    fig.text(.11, .035, 'PF50: three tidal rows, then signed arithmetic fusion of three channels. Same scale in all panels.\nZoom fixed from unshifted S1; the response is not labelled as confirmed fault contact. Overview retains block extrema.', fontsize=6, color='.3')
    save(fig, out, '01_signed_tsa_waveforms', qa)
    dump_csv(out/'waveform_display_data.csv', display)

    base_amp = data['none_rank0_amplitude']
    band = np.flatnonzero((freq >= 4000) & (freq <= 10000))
    spectral_center = float(freq[band[np.argmax(base_amp[band])]])
    spectrum_display = []
    fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 125/25.4))
    fig.subplots_adjust(left=.11, right=.98, top=.87, bottom=.15, hspace=.57, wspace=.38)
    for r, rank in enumerate([0, 1]):
        for c, ax in enumerate(axes[r]):
            for method, label, color in zip(METHODS, LABELS, COLORS):
                y = data[f'{method}_rank{rank}_amplitude']
                if c == 0:
                    ix = peak_preserving_indices(y, 0, np.searchsorted(freq, 12000), 2400)
                    x = freq[ix]/1000
                else:
                    ix = np.flatnonzero(abs(freq-spectral_center) <= 30)
                    x = freq[ix]
                ax.plot(x, y[ix], color=color, lw=.7, label=label)
                spectrum_display.extend(dict(panel=f'{r},{c}', method=method, rank=rank,
                                             frequency_hz=float(f), amplitude=float(v)) for f, v in zip(freq[ix], y[ix]))
            ax.set(xlabel='Equivalent frequency (kHz)' if c == 0 else 'Equivalent frequency (Hz)',
                   ylabel='Direct Fourier amplitude', xlim=(0, 12) if c == 0 else (spectral_center-30, spectral_center+30))
            mark(ax, chr(97+2*r+c), ('Direct TSA' if rank == 0 else 'Rank-1 SVD + TSA')+(' | spectrum' if c == 0 else ' | detail'))
    maximum = max(ax.get_ylim()[1] for ax in axes.flat)
    for ax in axes.flat:
        ax.set_ylim(0, maximum)
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.55, .985), ncol=3)
    fig.text(.11, .035, 'Direct FFT of signed fused TSA; 4-10 kHz input band, equivalent Fourier grid ~0.064 Hz.\nAll spectral axes start at zero. Zoom fixed from the unaligned spectrum; full-resolution arrays are retained.', fontsize=6, color='.3')
    save(fig, out, '02_signed_direct_spectra', qa)
    dump_csv(out/'spectrum_display_data.csv', spectrum_display)

    lags = read_csv(out/'lag_fits.csv')
    fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 85/25.4))
    fig.subplots_adjust(left=.105, right=.98, top=.82, bottom=.23, wspace=.35)
    for ax, channel in zip(axes, [2, 3]):
        for method, label, color in zip(METHODS[1:], LABELS[1:], COLORS[1:]):
            rows = [r for r in lags if r['condition'] == 'PF50' and r['method'] == method and int(r['channel']) == channel]
            values = [float(r['lag_ms']) for r in rows]
            if len(values) == 1:
                values *= 3
            ax.plot([1, 2, 3], values, color=color, marker='o', label=label)
        ax.axhline(0, color='.6', lw=.5)
        ax.set(xlabel='Full tidal row', ylabel='Signed waveform lag (ms)', xticks=[1, 2, 3])
        mark(ax, chr(95+channel), f'S{channel} relative to fixed S1')
    hh, ll = axes[0].get_legend_handles_labels()
    fig.legend(hh, ll, loc='upper center', bbox_to_anchor=(.55, .98), ncol=2)
    fig.text(.105, .035, 'One scalar time shift per selected span; Pearson maximum is an in-sample fitting criterion.\nNo individual-event registration or independent fault-contact timing is used.', fontsize=6, color='.3')
    save(fig, out, '03_signed_waveform_lags', qa)
    dump_json(out/'figure_qa.json', dict(backend='Python', figures=qa, visual_inspection='pending',
                                       waveform_zoom_center_s=center, spectrum_zoom_center_hz=spectral_center))


def report(out):
    rows = read_csv(out/'metrics.csv')
    def get(condition, method, rank):
        return next(r for r in rows if r['condition'] == condition and r['method'] == method and int(r['rank']) == rank)
    text = ['# 有符号振动波形对齐与同步平均预试验', '',
            '## 实际处理', '',
            '51,200 Hz、4-10 kHz 振动波形 → 波形 Pearson 时移估计 → 共同角域潮汐矩阵 → 每通道固定阶数 SVD/直接 TSA → 三通道有符号平均 → 直接傅里叶谱。匹配、重构、平均和频谱输入均未取包络。',
            '保持三组：未对齐、整段对齐、逐潮汐周期对齐。每个通道对齐到固定 S1；本轮没有额外加入故障事件定位或跨潮汐行配准。',
            '原稿 PSI 依赖模态包络，本轮按预试验约定暂不使用；秩 0 表示直接 TSA，秩 1/2 仅为固定降噪敏感性对照，不是 PSI-TPSVD。', '',
            '## PF50 结果', '',
            '| SVD 秩 | 方法 | 有符号融合 RMS | 通道相干能量比例 | 故障间隔梳状谱/邻域背景 dB | 啮合线能量比例 |',
            '|---:|---|---:|---:|---:|---:|']
    for rank in RANKS:
        for method, name in zip(METHODS, CN):
            r = get('PF50', method, rank)
            text.append(f'| {rank} | {name} | {float(r["signed_fused_rms"]):.5f} | {float(r["channel_coherent_energy_fraction"]):.4f} | {float(r["fault_spaced_comb_background_db"]):.3f} | {float(r["mesh_bin_energy_fraction"]):.4f} |')
    text += ['', '## 相对未对齐的增量', '', '| 状态 | 秩 | 方法 | 相干比例变化 | 梳状谱/背景变化 dB |', '|---|---:|---|---:|---:|']
    for condition in ['BL', 'PF50']:
        for rank in RANKS:
            base = get(condition, 'none', rank)
            for method, name in zip(METHODS[1:], CN[1:]):
                r = get(condition, method, rank)
                dc = float(r['channel_coherent_energy_fraction'])-float(base['channel_coherent_energy_fraction'])
                db = float(r['fault_spaced_comb_background_db'])-float(base['fault_spaced_comb_background_db'])
                text.append(f'| {condition} | {rank} | {name} | {dc:+.4f} | {db:+.3f} |')
    text += ['', '## 如何解释', '',
             '- 相干能量比例 = 三通道平均波形能量 / 三通道平均能量，范围 0-1；增大表示通道间更一致，不等于确认故障脉冲对齐。所有通道保留原极性，未为了增强结果翻转极性或补偿幅值。',
             '- 梳状谱指标在 4-10 kHz 中统计故障周期对应的高次谱线，不是在约 5.38 Hz 处测低频峰。目标为每 84 个 Fourier bin 的位置及左右各 1 bin，背景为同一位置左右 5-10 bin。排除啮合线附近 ±3 个行星架阶次并保留边界余量；使用 RMS 幅值比。',
             '- 该梳状谱指标是探索性周期结构指标，不是物理 SNR，也不能排除其他同周期来源。健康组同样计算，用于观察非故障特异的同步增强。',
             '- 波形减小可以来自噪声抵消，也可以来自故障被抵消；波形增大可以来自故障，也可以来自同步啮合。不能仅凭幅值选择优胜方法。',
             '- 按当前定义，啮合频率是故障频率的 31 倍，故障周期/潮汐周期平均均可能保留同步啮合。没有证据支持“除故障之外全部衰减”。',
             '- 相位与矩阵沿用当前齿数 21/31/84，几何一致性待核实。三个潮汐周期约 15.62 s/行，属于同一历史记录，不是三个独立实验。',
             '- 原始工作采样率为 51,200 Hz；角域行有 524,288 点，频率坐标是按平均行时长换算的等效频率。匹配窗口仍为约 ±8 ms，每个长片段只估计一个时移。',
             '- 本次从保存的带通振动波形重新计算时移、重采样和 SVD/TSA，并核对历史有符号分支数组；历史数据未修改。频域新结果不与以前的包络谱 dB 指标直接作优劣比较。', '',
             '完整数字见 metrics.csv；时移与相关曲线见 lag_fits.csv、lag_curves.csv；两个状态的完整三通道 TSA、融合波形和频谱保存在 signed_results.npz；图中显示点另有 CSV。', '']
    (out/'有符号波形预试验报告.md').write_text('\n'.join(text), encoding='utf-8')
    for rank in RANKS:
        print('PF50 rank', rank, [(m, get('PF50', m, rank)['channel_coherent_energy_fraction'],
                                   get('PF50', m, rank)['fault_spaced_comb_background_db']) for m in METHODS], flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, default=ROOT/'results'/'event_alignment_v1_20260906_053745')
    parser.add_argument('--prior', type=Path, default=ROOT/'results'/'pearson_tidal_vs_whole_20260906_092807')
    parser.add_argument('--plot-only', type=Path)
    args = parser.parse_args()
    if args.plot_only:
        out = args.plot_only.resolve()
    else:
        out = ROOT/'results'/('signed_waveform_pilot_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
        out.mkdir()
        print('OUTPUT', out, flush=True)
        compute(args.source.resolve(), args.prior.resolve(), out)
    report(out)
    plot(out)
    print('COMPLETE', out, flush=True)


if __name__ == '__main__':
    main()
