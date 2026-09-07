"""Arrival-focused figures and Chinese evidence report from frozen V2 outputs."""
from pathlib import Path
import argparse
import json
from plot_fault_event_alignment import np, plt, read_csv, pick, load_npz, mark, save, ROOT

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['svg.fonttype'] = 'none'
METHODS = ['none','waveform_xcorr','v1_constraint','envelope_xcorr','joint_envelope']
SHORT = ['Unaligned','Waveform\nXCorr','V1\nconstraint','Envelope\nXCorr','Joint\nenvelope']
ZH = dict(none='未对齐',waveform_xcorr='原始波形互相关',v1_constraint='上一版约束时移',
          envelope_xcorr='局部包络互相关',joint_envelope='新版联合包络匹配')
COLOR = dict(zip(METHODS,['#89919A','#8A91AB','#C5A58B','#648EA3','#B16D52']))


def avg(rows, field, **where):
    values=pick(rows,field,**where)
    return float(np.nanmean(values)) if np.isfinite(values).any() else float('nan')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('result',type=Path)
    root=parser.parse_args().result.resolve()
    manifest=json.loads((root/'run_manifest.json').read_text(encoding='utf-8'));assert manifest['complete']
    output=root/'figures';output.mkdir(exist_ok=True)
    metrics=read_csv(root/'metrics.csv');arrival=read_csv(root/'heldout_prediction.csv')
    ds=read_csv(root/'wrong_order_gain.csv');strata=read_csv(root/'envelope_strata.csv')
    spectra=read_csv(root/'spectra.csv');qa=[]

    fig,axes=plt.subplots(2,2,figsize=(183/25.4,145/25.4))
    fig.subplots_adjust(left=.13,right=.98,bottom=.15,top=.92,wspace=.43,hspace=.68)
    for k,method in enumerate(METHODS):
        values=pick(metrics,'cross_band_arrival_spread_ms',condition='PF50',method=method)
        axes[0,0].scatter(k+np.linspace(-.06,.06,len(values)),values,color=COLOR[method],s=18)
        axes[0,0].plot([k-.14,k+.14],[np.nanmean(values)]*2,color=COLOR[method],lw=1.5)
        rr=[r for r in strata if r['condition']=='PF50' and r['method']==method and r['hypothesis']=='target' and r['mask']=='all_five']
        pairs=sorted(set((r['block1'],r['block2']) for r in rr))
        vals=[np.mean([float(r['phase']) for r in rr if (r['block1'],r['block2'])==pair]) for pair in pairs]
        axes[0,1].scatter(k+np.linspace(-.06,.06,len(vals)),vals,color=COLOR[method],s=18)
        if vals:axes[0,1].plot([k-.14,k+.14],[np.mean(vals)]*2,color=COLOR[method],lw=1.5)
    for ax in axes[0]:ax.set(xticks=range(5),xticklabels=SHORT)
    axes[0,0].set(ylabel='Cross-band arrival spread (ms)',ylim=(0,None))
    axes[0,1].set_ylabel('Envelope consistency');axes[0,1].axhline(0,color='.7',lw=.6)
    mark(axes[0,0],'a','PF50: separate-band arrival proxy')
    mark(axes[0,1],'b','PF50: repeated packet concentration')
    for k,method in enumerate(['none','envelope_other_median','joint_other_prediction']):
        values=pick(arrival,'mae_ms',condition='PF50',method=method)
        color=['#89919A','#648EA3','#B16D52'][k]
        axes[1,0].scatter(k+np.linspace(-.07,.07,len(values)),values,color=color,s=15)
        axes[1,0].plot([k-.14,k+.14],[np.nanmean(values)]*2,color=color,lw=1.5)
    axes[1,0].set(xticks=range(3),xticklabels=['Nominal','Other-channel\nmedian','Joint\nprediction'],ylabel='Held-out onset proxy MAE (ms)',ylim=(0,None))
    mark(axes[1,0],'c','PF50: predict the third sensor')
    for condition,dx,color,marker in [('BL',-.09,'#648EA3','o'),('PF50',.09,'#B16D52','s')]:
        for k,method in enumerate(['envelope_xcorr','joint_envelope']):
            value=avg(ds,'delta_envelope_gain',condition=condition,method=method)
            axes[1,1].scatter(k+dx,value,color=color,marker=marker,s=25,label=condition if k==0 else None)
    axes[1,1].set(xticks=[0,1],xticklabels=['Envelope XCorr','Joint envelope'],ylabel='Target-minus-control gain')
    axes[1,1].axhline(0,color='.7',lw=.6);axes[1,1].legend(fontsize=6)
    mark(axes[1,1],'d','Matched wrong-order controls')
    fig.text(.13,.025,'a: three time blocks; b: three dependent block-pair means; c: three sensors × three blocks, common supported events.\nCross-band and held-out onsets are response proxies, not mechanical contact truth. One previously inspected recording per state.',fontsize=6,color='.3')
    save(fig,output,'01_arrival_results',qa)

    data=load_npz(root/'PF50'/'evaluated_envelopes.npz')
    bins=(np.mod(data['psi'],2*np.pi)/(np.pi/3)).astype(int)
    use=data['common_mask'] & (data['branch']==0) & (bins==0) & (data['t']>=27) & (data['t']<37)
    ids=np.flatnonzero(use)[:12]
    if not len(ids):raise ValueError('Fixed display stratum has no valid event; do not substitute a best-looking group.')
    fig,axes=plt.subplots(3,3,figsize=(183/25.4,151/25.4))
    fig.subplots_adjust(left=.12,right=.985,bottom=.15,top=.925,wspace=.36,hspace=.57)
    shown=['none','envelope_xcorr','joint_envelope'];u=data['u']*1000
    example_rows=[]
    for j in range(3):
        maximum=max(np.max(data[m][ids,j]) for m in shown)*1.05
        for k,method in enumerate(shown):
            ax=axes[j,k];yy=data[method][ids,j]
            ax.plot(u[::4],yy[:,::4].T,color=COLOR[method],alpha=.22,lw=.5)
            ax.plot(u[::4],yy.mean(axis=0)[::4],color=COLOR[method],lw=1.2)
            ax.axvline(0,color='.6',lw=.55,ls=':')
            ax.set(xlim=(-12,24),ylim=(0,maximum),xticks=[-10,0,10,20],xlabel='Local time (ms)')
            if k==0:ax.set_ylabel(f'S{j+1} envelope')
            mark(ax,chr(97+j*3+k),['Unaligned','Envelope XCorr','Joint envelope'][k] if j==0 else f'S{j+1}')
            example_rows += [dict(method=method,channel=j+1,time_ms=tt,mean_envelope=vv,event_count=len(ids))
                             for tt,vv in zip(u[::4],yy.mean(axis=0)[::4])]
    fig.text(.12,.025,f'PF50; fixed branch A, carrier angle 0–60°, 27–37 s; same {len(ids)} events. Thin: individual envelopes; thick: mean.\nWide windows retain the original packets. Same amplitude limits within each sensor; acquisition units, before SVD.',fontsize=6,color='.3')
    save(fig,output,'02_envelope_packets',qa)
    from run_event_arrival_v2 import save_csv
    save_csv(root/'display_source_data.csv',example_rows)

    fig,axes=plt.subplots(1,2,figsize=(183/25.4,104/25.4),gridspec_kw={'width_ratios':[1.5,1]})
    fig.subplots_adjust(left=.105,right=.985,bottom=.23,top=.88,wspace=.40)
    for k,method in enumerate(shown):
        x=pick(spectra,'carrier_order',condition='PF50',method=method,block=0)
        y=pick(spectra,'amplitude',condition='PF50',method=method,block=0)
        cut=x<=15
        axes[0].plot(x[cut],y[cut],color=COLOR[method],label=['Unaligned','Envelope XCorr','Joint envelope'][k],ls=['-','--','-'][k],lw=.85)
        vals=pick(metrics,'family_db',condition='PF50',method=method)
        axes[1].scatter(k+np.linspace(-.06,.06,len(vals)),vals,color=COLOR[method],s=19)
        axes[1].plot([k-.12,k+.12],[np.mean(vals)]*2,color=COLOR[method],lw=1.5)
    for k in range(1,6):axes[0].axvline(k*84/31,color='.7',lw=.55,ls=':',zorder=0)
    axes[0].set(xlim=(0,15),ylim=(0,None),xlabel='Carrier order',ylabel='Mean envelope amplitude')
    axes[0].legend(fontsize=5.8);mark(axes[0],'a','Continuous 27–37 s block')
    axes[1].set(xticks=range(3),xticklabels=['Unaligned','Envelope\nXCorr','Joint\nenvelope'],ylabel='Fault family / background (dB)')
    mark(axes[1],'b','PF50 diagnostic feature retention')
    fig.text(.105,.04,'Original time axis after bounded local time shifts; mean of individual channel envelopes, before event SVD.\nDotted guides: multiples of 84/31 carrier order. Three dots: three time blocks, not independent runs.',fontsize=6,color='.3')
    save(fig,output,'03_envelope_spectrum',qa)
    report(root,metrics,arrival,ds,data,manifest)
    (root/'figure_qa.json').write_text(json.dumps(dict(backend='Python',figures=qa,visual_inspection='pending'),indent=2),encoding='utf-8')
    print(f'OUTPUT={output}',flush=True)


def report(root,metrics,arrival,ds,data,manifest):
    lines=['# 故障脉冲到达时刻对齐：V2 真实数据试验','',
           '本轮目标是同一次冲击的事件时刻对齐，不要求内部高频振铃同相。算法使用振动信号；编码器未进入本轮匹配。',
           '沿用先前已查看的 BL/PF50 各一条记录、三通道 0/120/240°，开发/验证/测试按原划分。本轮为探索性内部留出，不是新的盲测或跨记录独立重复。','',
           '## 本轮结论','',
           '新版能够明显提高重复脉冲包络的集中程度，相比 V1 有改进；这一收益主要由局部包络事件匹配贡献，新增周期/空间约束尚未显示稳定的总体优势。',
           f'PF50 包络一致性：未对齐 {avg(metrics,"envelope_consistency",condition="PF50",method="none"):.3f}，包络互相关 {avg(metrics,"envelope_consistency",condition="PF50",method="envelope_xcorr"):.3f}，联合匹配 {avg(metrics,"envelope_consistency",condition="PF50",method="joint_envelope"):.3f}。',
           f'留出通道预测 MAE：简单其他通道基线 {avg(arrival,"mae_ms",condition="PF50",method="envelope_other_median"):.3f} ms，联合预测 {avg(arrival,"mae_ms",condition="PF50",method="joint_other_prediction"):.3f} ms，数值几乎相同。',
           '联合匹配的跨通道到达差比包络互相关小，但未优于未处理信号；因此不能把包络叠加图变尖直接写成“三通道同一故障事件已全部对齐”。目标与错误阶次的增益有正向差异，但相对局部包络互相关的额外收益很小，且没有独立重复证据。',
           '当前可支持的定位是：事件包络同步的可行版本已实现；物理约束带来的独立贡献及实测故障事件身份仍待验证。','',
           '## 实际实现','',
           '1. 沿用 4–10 kHz 带通响应和振动运动坐标；51200 Hz 支路取 Hilbert 幅值包络，Gaussian 平滑 σ=0.12 ms。',
           '2. 在校准段用局部背景以上 20% 包络上升沿作为到达标记；以其为零点，按分支/角度/通道构造单位范数事件模板。事件不足时回退同分支模板并记录。',
           '3. 局部包络归一化互相关搜索 ±8 ms，保留最多四个相隔至少 0.8 ms 的峰，最低相关系数 0.25，峰附近抛物线插值。',
           '4. 从校准段拟合零通道均值的相对周期路径；动态规划结合候选匹配代价、相对路径与公共事件偏差连续性，在候选之间选择。最终修正使用选中候选的时移，而不是全部替换成平滑模型值。',
           f'5. 三个约束强度仅按验证段留一通道误差及覆盖率选型，冻结强度 {manifest["selected_strength"]} 后计算测试；未根据本轮测试结果调参。','',
           '## PF50 结果','',
           '| 方法 | 包络一致性↑ | 跨频带通道到达差 ms↓ | 同频带上升沿到达差 ms↓ | 覆盖率 | 故障族/背景 dB |',
           '|---|---:|---:|---:|---:|---:|']
    for method in METHODS:
        lines.append(f'| {ZH[method]} | {avg(metrics,"envelope_consistency",condition="PF50",method=method):.4f} | '
                     f'{avg(metrics,"cross_band_arrival_spread_ms",condition="PF50",method=method):.3f} | '
                     f'{avg(metrics,"inband_arrival_spread_ms",condition="PF50",method=method):.3f} | '
                     f'{avg(metrics,"coverage",condition="PF50",method=method):.1%} | '
                     f'{avg(metrics,"family_db",condition="PF50",method=method):.2f} |')
    lines += ['',
        '包络一致性为同通道/同分支/同角度分层的跨块归一化包络内积，属于包络形状集中程度，不是故障识别准确率。此处宽窗 [-12,24) ms，与 V1 的 20 ms 窗数值不可直接跨报告比较。上表五方法一致性与到达差使用五方法交集事件；主三方法另有更大公共集合记录。',
        '跨频带到达差：对预先固定的 0.8–2.5 kHz 上升沿标记，减去各方法在 4–10 kHz 推断的时移，计算可用通道对之间的绝对时间差，再取块均值。不同频带有传播色散，故仅作旁证，不是真实机械接触误差。同频带上升沿虽未被测试段互相关直接优化，但仍共享信号证据，独立性弱于留一通道预测。',
        '覆盖率是候选事件被处理的比例，不是检测准确率；频谱为原时间轴上重构后的平均单通道包络，未串接短窗，未施加幅值逆补偿。','',
        '## 留一通道预测：相同评分事件','',
        '| 状态 | 未修正偏差 ms | 其他通道中值基线 ms | 联合预测偏差 ms | 联合预测自身覆盖率 | 公共评分覆盖率 |',
        '|---|---:|---:|---:|---:|---:|']
    for condition in ['BL','PF50']:
        vals=[avg(arrival,'mae_ms',condition=condition,method=m) for m in ['none','envelope_other_median','joint_other_prediction']]
        cov=avg(arrival,'own_coverage',condition=condition,method='joint_other_prediction')
        common=avg(arrival,'common_coverage',condition=condition,method='joint_other_prediction')
        lines.append(f'| {condition} | {vals[0]:.3f} | {vals[1]:.3f} | {vals[2]:.3f} | {cov:.1%} | {common:.1%} |')
    lines += ['',
        '中值基线：从校准段估计固定通道偏置，用其他两个通道的局部包络匹配时移预测第三通道。仅剩两个通道时中值等于均值。联合预测用冻结周期路径与其他通道候选序列，测试时不读取被留出通道的候选。参考为该通道预先独立固定的最强突出波包上升沿，仍是响应代理量，不是外部真值；标记也可能落在错误波包。',
        '表中平均为三个通道×三个块的 MAE 等权平均，非九次独立重复；两种预测方法在同一交集事件上评分，另列覆盖率。','',
        '## 错误阶次与健康控制','',
        '| 状态 | 方法 | 目标增益减错误阶次增益 | 匹配分层块对数 |','|---|---|---:|---:|']
    for row in ds:
        lines.append(f'| {row["condition"]} | {ZH[row["method"]]} | {float(row["delta_envelope_gain"]):.4f} | {row["matched_stratum_pairs"]} |')
    lines += ['',
        '错误阶次为目标事件阶次的 0.85/1.15 倍，各自用校准段生成事件、模板和路径，保持相同匹配自由度。上表使用主三方法公共集合，并在目标/两种错阶次中匹配分支、角度、通道和块对。108 个分层块对不是 108 次独立试验。健康记录也可有正常瞬态，此处没有建立故障分类器，因此不报告“健康检测准确率”。','',
        '## 已知时延控制测试','',
        '三个随机种子，已知事件计划，随机内部振铃相位与幅值、不同通道振铃频率及噪声。包络互相关和联合包络匹配的平均绝对时刻误差均约 0.029 ms，未对齐约 1.459 ms，测试事件覆盖 100%。本测试中二者相同，没有证明约束的额外优势。',
        '通过：留出通道测试数据被扰动时预测不变；校准路径冻结；全无证据时拒绝；相对路径零均值。这个简单单事件波包测试不含复杂故障源识别或多体动力学，不能替代实测故障身份验证。','',
        '## 交付与边界','',
        '- `metrics.csv`、`envelope_strata.csv`、`heldout_prediction.csv`、`wrong_order_gain.csv` 与 `spectra.csv` 为图表数值。',
        '- 各状态目录保存包络、模板、上升沿标记、候选时差、冻结路径、错误阶次结果及对齐包络数组。',
        '- `run_manifest.json` 保存实际参数及核心代码散列，`selection.json` 保存开发选型依据。',
        '- 主结果均在 SVD 前；这轮未把 V1 的 SVD 降噪重新运行作为新方法收益。',
        '- 原始波形互相关与 V1 结果来自上一轮，作为补充；本轮主要公平对照是共享模板与事件自由度的局部包络互相关。',
        '- 三行星/三传感器对称条件下仍有公共相位不可辨识性，不据本轮声称唯一分离了速度误差和真实传播时延。',
        '- 本轮评估的是事件对齐及诊断特征，尚无单个实测脉冲的外部故障身份真值。','']
    lines += ['## 后续最具体的检查','',
              '检查发生候选切换或通道到达差较大的事件：局部包络互相关是否匹配到了邻近正常啮合波包，或固定上升沿代理是否选错了最强波包。本版瞬态证据仅作为每事件/通道的可用性门槛，尚未对每一个互相关候选分别检验其异常能量证据；这应作为下一版明确的候选身份检查，而不是直接增加振铃同相约束。',
              '本轮没有在看到测试结果后启用该扩展，避免把事后修改与当前冻结试验混在一起。','']
    (root/'故障脉冲到达对齐_V2报告.md').write_text('\n'.join(lines),encoding='utf-8')


if __name__=='__main__':main()
