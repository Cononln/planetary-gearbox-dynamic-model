"""Python publication figures and an evidence-based Chinese pilot report."""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'runtime'/'phase_python'))
os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'runtime'/'matplotlib_feasibility'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'pdf.fonttype': 42, 'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .7,
                     'legend.frameon': False, 'legend.fontsize': 6, 'lines.linewidth': 1.0})

COL = {'none': '#7A828C', 'local_xcorr': '#7A9CAF', 'constrained': '#AA6D58',
       'BL': '#7895A8', 'PF50': '#AA6D58'}
LABEL = {'none': 'Unaligned', 'local_xcorr': 'Local XCorr', 'constrained': 'Event-constrained'}
METHODS = list(LABEL)


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def pick(rows, field, **where):
    return np.asarray([float(r[field]) for r in rows if all(r[k] == str(v) for k, v in where.items())])


def load_npz(path):
    with np.load(path) as z:
        return {key: z[key] for key in z.files}


def mark(ax, letter, title):
    ax.text(-.16, 1.05, letter, transform=ax.transAxes, fontweight='bold', fontsize=9)
    ax.set_title(title, loc='left', pad=7)


def save(fig, out, name, qa):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    clipped = []
    off_view_ticks = set()
    for ax in fig.axes:
        for axis in [ax.xaxis, ax.yaxis]:
            lo, hi = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if not lo-1e-10 <= tick.get_loc() <= hi+1e-10:
                    off_view_ticks.update([tick.label1, tick.label2])
    for item in fig.findobj(matplotlib.text.Text):
        if item in off_view_ticks or not item.get_visible() or not item.get_text().strip():
            continue
        box = item.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width+1 or box.y1 > fig.bbox.height+1:
            clipped.append(item.get_text())
    for ext in ['png', 'svg', 'pdf', 'tiff']:
        args = {'dpi': 220 if ext == 'png' else 600}
        if ext == 'tiff':
            args['pil_kwargs'] = {'compression': 'tiff_lzw'}
        fig.savefig(out/f'{name}.{ext}', **args)
    xml = ET.parse(out/f'{name}.svg')
    texts = len(xml.findall('.//{http://www.w3.org/2000/svg}text'))
    assert texts > 10
    qa.append(dict(figure=name, svg_text_elements=texts, width_mm=round(fig.get_figwidth()*25.4, 1),
                   height_mm=round(fig.get_figheight()*25.4, 1), formats=['svg', 'pdf', 'png', 'tiff'],
                   text_outside_canvas=clipped))
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('result', type=Path)
    args = parser.parse_args(); root = args.result.resolve()
    status = json.loads((root/'evaluation_status.json').read_text(encoding='utf-8'))
    assert status['complete']
    output = root/'figures'; output.mkdir(exist_ok=True)
    cfg = json.loads((root/'runtime_config.json').read_text(encoding='utf-8'))
    selected = json.loads((root/'selection.json').read_text(encoding='utf-8'))
    metrics = read_csv(root/'metrics.csv'); phases = read_csv(root/'phase_strata.csv')
    spectra = read_csv(root/'spectra.csv'); arrival = read_csv(root/'arrival_prediction.csv')
    delta = read_csv(root/'delta_S.csv'); qa = []
    data = load_npz(root/'PF50'/'evaluated_events.npz')
    folder = root/'PF50'/f'band_{selected["band_index"]}'
    motion = load_npz(root/'PF50'/'motion.npz')
    evidence = np.load(folder/'evidence.npy', mmap_mode='r')
    wave = np.load(folder/'waveform.npy', mmap_mode='r')
    ev = load_npz(folder/'events.npz'); ca = load_npz(folder/'candidates.npz')
    fhi, flo = cfg['working_rates_hz']['impulse'], cfg['working_rates_hz']['motion']
    # Fixed display interval, not selected for largest alignment gain.
    span = (28.0, 28.5)
    fig, axes = plt.subplots(3, 2, figsize=(183/25.4, 155/25.4))
    fig.subplots_adjust(left=.13, right=.985, bottom=.125, top=.93, wspace=.38, hspace=.67)
    for j in range(3):
        a, b = [int(s*fhi) for s in span]
        tt = np.arange(a, b, 4)/fhi
        axes[j, 0].plot(tt, wave[a:b:4, j], color='#536B82', lw=.55)
        a, b = [int(s*flo) for s in span]
        axes[j, 1].plot(np.arange(a, b)/flo, evidence[a:b, j], color='#7D8A96', lw=.65)
        axes[j, 1].axhline(3, color=COL['constrained'], ls=':', lw=.7)
        for n in np.flatnonzero((ev['t'] >= span[0]) & (ev['t'] < span[1])):
            for ax in axes[j]:
                ax.axvline(ev['t'][n], color=['#AA6D58', '#587F9E'][ev['branch'][n]], alpha=.65, lw=.7, ls='--')
            if np.isfinite(ca['offsets'][n, j, 0]):
                tx = ev['t'][n]+ca['offsets'][n, j, 0]
                axes[j, 1].scatter(tx, np.interp(tx, motion['time'], evidence[:, j]), s=13, color=COL['constrained'], zorder=5)
        axes[j, 0].set(ylabel='Acquisition amplitude', xlim=span, xlabel='Time (s)')
        axes[j, 1].set(ylabel='Transient contrast D', xlim=span, xlabel='Time (s)')
        mark(axes[j, 0], chr(97+2*j), f'S{j+1}: band response')
        mark(axes[j, 1], chr(98+2*j), f'S{j+1}: event candidates')
    fig.text(.13, .02, 'PF50; fixed 28.0–28.5 s window. Dashed lines: inferred A/B event schedules; dots: strongest precomputed candidates.\nCandidate peaks are not independently verified fault contacts. Band: %g–%g Hz.' % tuple(selected['common_band_hz']), fontsize=6, color='.3')
    save(fig, output, '01_event_detection', qa)

    bins = (np.mod(data['psi'], 2*np.pi)/(np.pi/3)).astype(int)
    use = data['common_mask'] & (data['branch'] == 0) & (bins == 0) & (data['t'] >= 27) & (data['t'] < 37)
    ids = np.flatnonzero(use)[:12]
    if not len(ids):
        raise ValueError('Predeclared event display stratum has no common-mask events; do not cherry-pick a replacement.')
    uu = data['local_u_s']*1000
    fig, axes = plt.subplots(3, 3, figsize=(183/25.4, 150/25.4))
    fig.subplots_adjust(left=.115, right=.985, bottom=.145, top=.93, wspace=.36, hspace=.52)
    waveform_source = []
    for j in range(3):
        limit = max(np.abs(data[m if m != 'none' else 'raw'][ids, j]).max() for m in METHODS)*1.05
        for k, m in enumerate(METHODS):
            ax = axes[j, k]; key = 'raw' if m == 'none' else m
            yy = data[key][ids, j]
            ax.plot(uu, yy.T, color=COL[m], alpha=.17, lw=.5)
            ax.plot(uu, yy.mean(axis=0), color=COL[m], lw=1.1)
            ax.set(xlim=(-4, 16), ylim=(-limit, limit), xlabel='Local time (ms)')
            if k == 0:
                ax.set_ylabel(f'S{j+1} amplitude')
            mark(ax, chr(97+j*3+k), LABEL[m] if j == 0 else f'S{j+1}')
            for t0, value in zip(uu[::4], yy.mean(axis=0)[::4]):
                waveform_source.append(dict(method=m, channel=j+1, local_time_ms=t0, mean_amplitude=value,
                                             events=len(ids), branch='A', angle_bin_deg='0-60', block='27-37'))
    fig.text(.115, .025, f'PF50, branch A, carrier angle 0–60°, 27–37 s; same {len(ids)} events in every column. Thin: events; thick: mean.\nSigned band waveforms are shown before final SVD. Identical per-channel amplitude limits; no per-method normalization.', fontsize=6, color='.3')
    save(fig, output, '02_pulse_alignment', qa)
    with (root/'display_waveform_means.csv').open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=list(waveform_source[0])); w.writeheader(); w.writerows(waveform_source)

    fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 105/25.4))
    fig.subplots_adjust(left=.10, right=.91, bottom=.22, top=.88, wspace=.48)
    a = data['branch'] == 0
    angle = np.rad2deg(np.mod(data['psi'][a], 2*np.pi)); sort = np.argsort(angle)
    for j, color in enumerate(['#587E9A', '#AF775E', '#7D8F79']):
        axes[0].plot(angle[sort], data['path_s'][a, j][sort]*1000, color=color, label=f'S{j+1}')
    axes[0].axhline(0, color='.7', lw=.6)
    axes[0].set(xlabel='Relative carrier angle (deg)', ylabel='Effective relative delay (ms)', xlim=(0, 360), ylim=(-4.2, 4.2))
    axes[0].legend(ncol=3, loc='upper right'); mark(axes[0], 'a', 'Frozen periodic delay model: branch A')
    test = np.zeros(len(data['t']), bool)
    for lo, hi in cfg['split_seconds']['test']:
        test |= (data['t'] >= lo) & (data['t'] < hi)
    zz = data['constrained_shift'][test]*1000
    zz = np.ma.array(zz, mask=np.broadcast_to(~data['constrained_mask'][test, None], zz.shape))
    im = axes[1].imshow(zz, aspect='auto', origin='lower', cmap='RdBu_r', vmin=-8, vmax=8, interpolation='nearest')
    axes[1].set(xticks=[0, 1, 2], xticklabels=['S1', 'S2', 'S3'], ylabel='Scheduled event index (test blocks)')
    cb = fig.colorbar(im, ax=axes[1], fraction=.045, pad=.04); cb.set_label('Total correction (ms)')
    mark(axes[1], 'b', 'Event corrections; blank = rejected')
    fig.text(.10, .045, 'Smooth model curves are fitted operational parameters, not measured structural transfer-path truth.\nThe total correction includes a common event residual. The reference is the zero-channel-mean delay.', fontsize=6, color='.3')
    save(fig, output, '03_delay_fields', qa)

    fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 140/25.4))
    fig.subplots_adjust(left=.115, right=.985, bottom=.13, top=.93, wspace=.4, hspace=.63)
    for m in METHODS:
        raw_key = 'raw' if m == 'none' else m
        # Show both single-channel and signed-fusion outcomes, not only an envelope.
        yy = data[m+'_svd'][ids]
        axes[0, 0].plot(uu, yy[:, 0].mean(axis=0), color=COL[m], label=LABEL[m])
        axes[0, 1].plot(uu, yy.mean(axis=(0, 1)), color=COL[m], label=LABEL[m])
        order = pick(spectra, 'carrier_order', condition='PF50', method=m, block=0)
        amp = pick(spectra, 'amplitude', condition='PF50', method=m, block=0)
        mask = (order >= 0) & (order <= 15)
        axes[1, 0].plot(order[mask], amp[mask], color=COL[m], lw=.9)
        value = pick(metrics, 'fp_family_db', condition='PF50', method=m)
        k = METHODS.index(m)
        axes[1, 1].scatter(k+np.linspace(-.07, .07, 3), value, color=COL[m], s=18)
        axes[1, 1].plot([k-.13, k+.13], [value.mean()]*2, color=COL[m], lw=1.5)
    for ax in axes[0]:
        ax.set(xlabel='Local time (ms)', ylabel='TSA amplitude', xlim=(-4, 16))
    axes[0, 0].legend(); mark(axes[0, 0], 'a', 'S1 TSA: fixed-rank local SVD')
    mark(axes[0, 1], 'b', 'Signed three-channel TSA after SVD')
    for k in range(1, 6):
        axes[1, 0].axvline(k*84/31, color='.65', lw=.55, ls=':', zorder=0)
    axes[1, 0].set(xlabel='Carrier order', ylabel='Mean envelope amplitude', xlim=(0, 15), ylim=(0, None))
    mark(axes[1, 0], 'c', 'Continuous 10 s block: envelope spectrum')
    axes[1, 1].set(xticks=range(3), xticklabels=['Unaligned', 'Local\nXCorr', 'Event\nconstraint'], ylabel='Fault family / background (dB)')
    mark(axes[1, 1], 'd', 'PF50 fault-family retention')
    fig.text(.115, .025, 'TSA: fixed event stratum shown in Fig. 2. Spectrum: original time axis, 27–37 s, mean of individual envelopes, before event SVD.\nDots: three time blocks from one PF50 recording; ticks: arithmetic mean. Dotted spectrum lines: multiples of 84/31 carrier order.', fontsize=6, color='.3')
    save(fig, output, '04_tsa_and_spectrum', qa)

    fig, axes = plt.subplots(2, 2, figsize=(183/25.4, 136/25.4))
    fig.subplots_adjust(left=.12, right=.985, bottom=.13, top=.93, wspace=.40, hspace=.64)
    for condition, offset in [('BL', -.11), ('PF50', .11)]:
        for k, m in enumerate(METHODS):
            cov = pick(metrics, 'coverage', condition=condition, method=m)
            axes[0, 0].scatter(k+offset+np.linspace(-.025, .025, 3), cov, color=COL[condition], s=15,
                               marker='o' if condition == 'BL' else 's', label=condition if k == 0 else None)
            r = [p for p in phases if p['condition'] == condition and p['method'] == m and p['hypothesis'] == 'target']
            pairs = sorted(set((p['block1'], p['block2']) for p in r))
            value = [np.mean([float(p['phase']) for p in r if (p['block1'], p['block2']) == pair]) for pair in pairs]
            axes[0, 1].scatter(k+offset+np.linspace(-.025, .025, len(value)), value, color=COL[condition], s=15,
                               marker='o' if condition == 'BL' else 's')
        for k, m in enumerate(['local_xcorr', 'constrained']):
            ds = pick(delta, 'delta_S', condition=condition, method=m)
            if len(ds) and np.isfinite(ds[0]):
                axes[1, 0].scatter(k+offset, ds[0], color=COL[condition], s=28, marker='o' if condition == 'BL' else 's')
            else:
                axes[1, 0].text(k+offset, 0, 'NA', ha='center', fontsize=6, color=COL[condition])
    for ax in axes[0]:
        ax.set(xticks=range(3), xticklabels=['Unaligned', 'Local\nXCorr', 'Event\nconstraint'])
    axes[0, 0].set(ylabel='Accepted candidate coverage', ylim=(0, 1.03)); axes[0, 0].legend()
    axes[0, 1].set_ylabel('Cross-block signed phase score')
    axes[0, 1].axhline(0, color='.7', lw=.6)
    mark(axes[0, 0], 'a', 'Coverage is reported with alignment')
    mark(axes[0, 1], 'b', 'Target-event waveform consistency')
    axes[1, 0].axhline(0, color='.65', lw=.65)
    axes[1, 0].set(xticks=[0, 1], xticklabels=['Local XCorr', 'Event constraint'], ylabel='Event-specific gain, delta S')
    mark(axes[1, 0], 'c', 'Matched wrong-order controls')
    alternatives = ['constrained', 'independent', 'no_path_smoothing', 'heldout_constrained']
    for k, m in enumerate(alternatives):
        val = pick(metrics, 'common_mask_phase', condition='PF50', method=m)
        if len(val):
            axes[1, 1].scatter(k, val[0], color=COL['PF50'], s=25)
    axes[1, 1].set(xticks=range(4), xticklabels=['Full', 'Free\nchannels', 'Weak\nsmoothing', 'Held-out\nchannel'], ylabel='PF50 signed phase score')
    axes[1, 1].axhline(0, color='.7', lw=.6); mark(axes[1, 1], 'd', 'Model variants and prediction check')
    fig.text(.12, .025, 'Coverage dots: three time blocks. Phase dots: three block-pair averages over compatible event strata.\nNo independent-repeat significance is inferred. Variant scores use supported intersections; differing coverage must be considered.', fontsize=6, color='.3')
    save(fig, output, '05_validation_and_controls', qa)

    fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 96/25.4))
    fig.subplots_adjust(left=.11, right=.985, bottom=.23, top=.88, wspace=.35)
    for ax, condition in zip(axes, ['BL', 'PF50']):
        for j in range(1, 4):
            for field, dx, color, label in [('uncorrected_proxy_mae_ms', -.11, '#87939E', 'Nominal event'),
                                           ('predicted_proxy_mae_ms', .11, COL[condition], 'Other-channel prediction')]:
                value = pick(arrival, field, condition=condition, channel=j)
                ax.scatter(j+dx+np.linspace(-.025, .025, len(value)), value, s=17, color=color, label=label if j == 1 else None)
                if np.isfinite(value).any():
                    ax.plot([j+dx-.065, j+dx+.065], [np.nanmean(value)]*2, color=color, lw=1.5)
        ax.set(xticks=[1, 2, 3], xticklabels=['S1 held out', 'S2 held out', 'S3 held out'],
               ylabel='Arrival proxy MAE (ms)', ylim=(0, 5))
        ax.legend(loc='lower right', fontsize=5.7)
        mark(ax, 'a' if condition == 'BL' else 'b', condition)
    fig.text(.11, .055, 'Each held-out channel is predicted from the other two and a frozen calibration path. Three dots = three time blocks.\nReference: strongest precomputed candidate in the held-out channel, not a measured physical fault-contact ground truth.', fontsize=6, color='.3')
    save(fig, output, '06_heldout_arrival', qa)
    diagnostic_file = root/'envelope_signed_strata.csv'
    if diagnostic_file.exists():
        diagnostic = read_csv(diagnostic_file)
        fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 103/25.4), sharey=True)
        fig.subplots_adjust(left=.11, right=.985, bottom=.24, top=.80, wspace=.25)
        all_values = []
        for ax, condition in zip(axes, ['BL', 'PF50']):
            for representation, dx, color, marker in [('signed', -.11, '#536B82', 'o'),
                                                       ('envelope', .11, '#AA6D58', 's')]:
                for k, method in enumerate(METHODS):
                    rr = [row for row in diagnostic if row['condition'] == condition and
                          row['method'] == method and row['representation'] == representation]
                    pairs = sorted(set((r['block1'], r['block2']) for r in rr))
                    values = [np.mean([float(r['phase']) for r in rr if (r['block1'], r['block2']) == pair])
                              for pair in pairs]
                    all_values.extend(values)
                    ax.scatter(k+dx+np.linspace(-.025, .025, len(values)), values, s=18, color=color, marker=marker,
                               label=representation.title() if k == 0 else None)
                    ax.plot([k+dx-.07, k+dx+.07], [np.mean(values)]*2, color=color, lw=1.5)
            ax.axhline(0, color='.7', lw=.65)
            ax.set(xticks=range(3), xticklabels=['Unaligned', 'Local\nXCorr', 'Event\nconstraint'])
            mark(ax, 'a' if condition == 'BL' else 'b', condition)
        axes[0].set_ylabel('Cross-block normalized coherence')
        axes[0].set_ylim(min(-.03, min(all_values)-.03), max(.45, max(all_values)+.03))
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.55, .95), ncol=2)
        fig.text(.11, .04, 'Post-hoc diagnostic; same accepted events and frozen time shifts, before SVD. Envelope: continuous-record Hilbert magnitude,\n0.527 ms smoothing. Dots: three dependent block-pair means; ticks: overall mean. No new fault-specificity control is implied.', fontsize=6, color='.3')
        save(fig, output, '07_envelope_vs_signed', qa)
    save_report(root, metrics, delta, arrival, selected, qa)
    print(f'FIGURES={output}', flush=True)


def save_report(root, metrics, delta, arrival, selected, qa):
    lines = ['# 故障事件约束时移：首轮真实数据实现报告', '',
             '已运行两条真实记录上的事件检测、约束时移、局部互相关、留一通道预测、局部 SVD/TSA 和错误阶次控制。',
             '本轮为同记录开发/验证/后续时间块的探索性内部留出；这些记录此前已被查看。每种状态只有一条记录，不给独立重复显著性结论。', '',
             f'开发段冻结的共同冲击频带：{selected["common_band_hz"][0]}–{selected["common_band_hz"][1]} Hz；冲击支路 51200 Hz，运动支路 5120 Hz。',
             '三测点为 0/120/240°；编码器未参与算法推断。当前实现只校正事件时移，额外振铃相位补偿关闭。', '',
             '## 留出结果', '',
             '| 状态 | 方法 | 测试事件覆盖率 | 同集合相干分数 Rφ | 故障频率族/背景 dB |',
             '|---|---|---:|---:|---:|']
    for condition in ['BL', 'PF50']:
        for method in METHODS:
            cov = pick(metrics, 'coverage', condition=condition, method=method).mean()
            phase = pick(metrics, 'common_mask_phase', condition=condition, method=method)[0]
            family = pick(metrics, 'fp_family_db', condition=condition, method=method).mean()
            lines.append(f'| {condition} | {LABEL[method]} | {cov:.1%} | {phase:.5f} | {family:.2f} |')
    lines += ['', '覆盖率中的“事件”是当前故障假设下的候选接触位置；健康记录也能有候选，不能把这一比例当检测准确率。',
              '相干分数在最终 SVD 前计算。Unaligned/Local XCorr/Constrained 共用同一目标事件掩码；频谱来自原时间轴上的平均单通道包络，未将短事件串接成新频谱。', '',
              'Rφ 统计同一通道、同事件分支、同行星架角度层内的跨块波形一致性，再跨层平均；它不是三通道相位差消失的直接测量。三通道有符号融合另作示例展示，不能替代跨通道相位验证。', '',
              '## 事件特异性相干增益 ΔS', '', '| 状态 | 方法 | ΔS | 匹配分层块对数 |', '|---|---|---:|---:|']
    for row in delta:
        v = float(row['delta_S'])
        value = f'{v:.5f}' if np.isfinite(v) else '未定义（匹配控制不足）'
        lines.append(f'| {row["condition"]} | {row["method"]} | {value} | {row["matched_stratum_pairs"]} |')
    pf_raw = pick(metrics, 'common_mask_phase', condition='PF50', method='none')[0]
    pf_new = pick(metrics, 'common_mask_phase', condition='PF50', method='constrained')[0]
    pf_local = pick(metrics, 'common_mask_phase', condition='PF50', method='local_xcorr')[0]
    ds_new = pick(delta, 'delta_S', condition='PF50', method='constrained')[0]
    ds_local = pick(delta, 'delta_S', condition='PF50', method='local_xcorr')[0]
    lines += ['', '## 直接判断', '',
              '当前约束时移版本没有通过故障脉冲波形对齐的可行性验证，不能写成已优于局部互相关。',
              f'PF50 原幅值归一化评分副本的相干变化：未对齐 {pf_raw:.5f}，约束时移 {pf_new:.5f}，局部互相关 {pf_local:.5f}。',
              f'错误阶次控制后的 ΔS：约束时移 {ds_new:.5f}，局部互相关 {ds_local:.5f}。',
              '这些数值需要联合解释：局部相关可能在健康/错误周期窗口上也产生较高波形一致性；候选时刻更接近也不保证高频有符号振铃同相。', '',
              '## 留一通道预测', '', '| 状态 | 通道 | 未校正候选偏差 ms | 冻结模型预测偏差 ms | 预测评分覆盖率 |', '|---|---|---:|---:|---:|']
    for condition in ['BL', 'PF50']:
        for j in [1, 2, 3]:
            before = pick(arrival, 'uncorrected_proxy_mae_ms', condition=condition, channel=j).mean()
            after = pick(arrival, 'predicted_proxy_mae_ms', condition=condition, channel=j).mean()
            cov = pick(arrival, 'coverage', condition=condition, channel=j).mean()
            lines.append(f'| {condition} | S{j} | {before:.3f} | {after:.3f} | {cov:.1%} |')
    lines += ['', '此处参考是被留出通道中预先固定的最强瞬态候选，不是外部测得的真实故障到达时刻。']
    diagnostic_file = root/'envelope_signed_diagnostic.csv'
    if diagnostic_file.exists():
        diagnostic = read_csv(diagnostic_file)
        lines += ['', '## 追加定位诊断：包络与有符号波形', '',
                  '这是本轮结果出来后的探索性诊断：复用已冻结的时移和相同事件掩码，不重新拟合，不调整事件或参数。',
                  '在完整带通信号上取 Hilbert 幅值包络并作 27 点（0.527 ms）居中平滑，再按相同的时移取窗，避免短窗 Hilbert 边界制造共同形状。两种表示均去均值、加窗、逐事件单位范数化，仅评分副本如此处理。', '',
                  '| 状态 | 方法 | 包络一致性 | 有符号波形一致性 |', '|---|---|---:|---:|']
        for condition in ['BL', 'PF50']:
            for method in METHODS:
                envelope = pick(diagnostic, 'coherence', condition=condition, method=method, representation='envelope')[0]
                signed = pick(diagnostic, 'coherence', condition=condition, method=method, representation='signed')[0]
                lines.append(f'| {condition} | {LABEL[method]} | {envelope:.5f} | {signed:.5f} |')
        lines += ['',
                  'PF50 上，约束时移改善了包络集中程度，但没有改善高频振铃的有符号波形一致性；局部互相关在这两项上都更高。BL 上约束后的包络一致性反而下降。因此本轮只能支持部分事件包络定位的可行性，不能支持完整故障脉冲相位对齐已完成。',
                  '包络分数没有新跑错误阶次对照，不能仅据健康/故障差异宣称故障特异性。这里 108 个分层块对来自每状态一条记录、三个时间块，非 108 次独立试验。', '',
                  '### 下一轮应验证的具体修改', '',
                  '1. 保留本版时移作为粗定位，在校准/验证段检查剩余偏差能否用单一精细时移解释；进入 51200 Hz 波形支路估计，输出微秒级残余时移及置信度。当前 4–10 kHz 的振荡周期仅为 100–250 μs，毫秒级候选时刻改善本身不足以保证同相；现有 5120 Hz 候选支路虽然有抛物线插值，也没有独立的高频相位观测。',
                  '2. 用校准段的同分支、同角度事件构造参考，以多频点相位随频率的斜率检查纯时移模型，保留周期/空间约束；在未用于拟合的频点或时间块检验残差。带宽不足、相干低或多峰歧义大时拒绝，不把每事件的任意相位旋转当作成功。',
                  '3. 只有单时移模型仍有可重复、可预测的结构化相位残差时，才考虑低自由度频率相关相位模型；若残差不稳定，先修正候选事件对应关系。输出三个独立判据：包络时刻、波形相干、目标事件相对于错误阶次的增益。',
                  '4. 扩展的全部选择只使用开发数据，保留 BL/错误阶次/局部互相关对照，在新的未参与选择的记录或时间块上确认。当前已查看的测试结果不能再次称为盲测。',
                  '以上是下一版的进入条件，本次未启用这些扩展，也未重新搜索参数以获得正结果。']
    lines += ['', '## 代码验证与适用边界', '',
              '- 已知候选时刻控制测试：时移 RMSE 从 1.184 ms 降到 0.082 ms；该测试验证约束估计，不验证从实测信号中发现故障事件。',
              '- 留出通道数据扰动不改变对该通道的预测；路径系数在测试阶段冻结；零证据拒绝；分数延时插值与零移位重构测试通过。',
              '- 初次求解的部分健康记录触碰迭代上限；原始失败结果保留。增加求解迭代上限后，重新进行开发选型，再运行本轮测试。',
              '- 局部事件 SVD 与完整 31 圈行星架潮汐矩阵不同，不能据此增加完整潮汐重复数量。',
              '- 全局互相关是单个行星相对转动周期的整数采样移位；主对照局部互相关具有亚采样估计。全局对照仅作补充。',
              '- 此处局部互相关为通用操作性基线，不等于已完整复现两篇参考论文的全部算法。',
              '- 优化历史中的 loss 为给定候选分配后的回归损失，不代表已证明完整离散/连续联合目标单调下降。',
              '- 时移可将原来位于窗口边缘或窗口外的脉冲移入核心窗口，图中峰值变大不等于幅值补偿或传播能量恢复。',
              '- 暂未将单个有效时延场解释成真实结构传播时延，也未声称唯一分离转速与公共路径相位。',
              '- 未使用机械零位或齿号；A/B 仅为推断出的事件序列，不能直接赋予已验证的物理接触标签。',
              '- 负对照实现了相邻啮合背景构造、BL 和两种错误事件阶次。相邻啮合背景未另外生成一个完整对齐重跑的评分控制组；ΔS 使用匹配错误阶次组。',
              '- 本版 SVD 为固定秩 3 的校准基投影；可用组不足时旁路并记录。尚未启用振铃相位扩展和跨记录独立确认。', '',
              '## 输出文件', '',
              '- `runtime_config.json`：实际参数、源缓存标识与代码 SHA256。',
              '- `development_candidates.csv` / `selection.json`：开发阶段选择依据。',
              '- `metrics.csv`、`phase_strata.csv`、`delta_S.csv`、`arrival_prediction.csv`、`coverage.csv`：原始数值。',
              '- `envelope_signed_diagnostic.csv` / `.json`、`envelope_signed_strata.csv`：追加包络/波形同集合诊断及方法说明。',
              '- 各状态 `evaluated_events.npz`：原始/校正/降噪事件、时移、路径和掩码；事件检测与原信号带通缓存另存于该状态目录。',
              '- `figures/`：六组主结果图及一组包络/波形诊断图（如已计算），含 SVG/PDF/PNG/TIFF；布局与数据说明见每图页脚。', '']
    (root/'真实数据实现报告.md').write_text('\n'.join(lines), encoding='utf-8')
    (root/'figure_qa.json').write_text(json.dumps(dict(backend='Python', figures=qa,
         data='Measured cached event arrays and exported CSV; no fabricated waveform measurements',
         visual_inspection='Pending model review of generated previews'), indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
