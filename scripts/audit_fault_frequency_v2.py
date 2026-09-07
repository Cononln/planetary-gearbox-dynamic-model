"""Post-hoc Hz spectrum and exact-order fault-harmonic audit of frozen V2 shifts."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0,str(ROOT/'src_py'))
import numpy as np
from scipy import signal
from event_alignment import reconstruct_block
from run_event_arrival_v2 import save_csv,save_json,load_npz
from plot_fault_event_alignment import plt,read_csv,pick,mark,save

plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans','Liberation Sans']
plt.rcParams['svg.fonttype']='none'
METHODS=['none','envelope_xcorr','joint_envelope']
LABEL=['Unaligned','Envelope XCorr','Joint envelope']
COLORS=['#858D96','#648EA3','#B16D52']


def coefficients(env,fs,start,motion,zr=84,zp=31):
    t=start+np.arange(len(env))/fs
    phase=np.interp(t,motion['time'],motion['phase'])/(2*np.pi*zr)
    phase-=phase[0]
    count=max(256,int(phase[-1]*4096))
    theta=np.linspace(0,phase[-1],count,endpoint=False)
    yy=np.column_stack([np.interp(theta,phase,env[:,j]) for j in range(3)])
    yy-=yy.mean(axis=0)
    hann=np.hanning(count)
    orders=np.arange(1,9)*zr/zp
    coeff=2*np.exp(-2j*np.pi*orders[:,None]*theta)@(yy*hann[:,None])/hann.sum()
    return coeff


def compute(root,out):
    manifest=json.loads((root/'run_manifest.json').read_text(encoding='utf-8'));assert manifest['complete']
    source=Path(manifest['source'])
    cfg=json.loads((source/'runtime_config.json').read_text(encoding='utf-8'))
    fs=cfg['working_rates_hz']['impulse'];zr=cfg['gear_teeth']['ring'];zp=cfg['gear_teeth']['planet']
    frequency,lines,curves,fusion=[],[],[],[]
    for condition in ['BL','PF50']:
        data=load_npz(root/condition/'evaluated_envelopes.npz')
        events={key:data[key] for key in ['t','psi','branch']}
        motion=load_npz(source/condition/'motion.npz')
        wave=np.load(source/condition/'band_2'/'waveform.npy',mmap_mode='r')
        for block,span in enumerate(manifest['split_seconds']['test']):
            fm=float(np.diff(np.interp(span,motion['time'],motion['phase']))[0]/(2*np.pi*(span[1]-span[0])))
            fc=fm/zr;fp=fm/zp
            frequency.append(dict(condition=condition,block=block,start_s=span[0],end_s=span[1],
                                  mesh_hz=fm,carrier_hz=fc,planet_same_side_hz=fp,two_contacts_hz=2*fp,
                                  time_fft_bin_hz=1/(span[1]-span[0])))
            for method in METHODS:
                print(f'{condition} block {block}: {method}',flush=True)
                y,derivative,rejections=reconstruct_block(wave,fs,events,data[method+'_shift'],data[method+'_mask'],span)
                assert derivative>=.199
                env=np.abs(signal.hilbert(y,axis=0))
                hann=np.hanning(len(env));centered=env-env.mean(axis=0)
                fft=2*np.fft.rfft(centered*hann[:,None],axis=0)/hann.sum()
                f=np.fft.rfftfreq(len(env),1/fs)
                average=np.abs(fft.mean(axis=1))
                keep=(f<=60)
                curves += [dict(condition=condition,block=block,method=method,frequency_hz=x,amplitude=a)
                           for x,a in zip(f[keep],average[keep])]
                exact=coefficients(env,fs,span[0],motion,zr,zp)
                for h in range(1,9):
                    use=np.flatnonzero(np.abs(f-h*fp)<=.18*fc)
                    ix=use[np.argmax(average[use])]
                    c=exact[h-1];incoherent=np.mean(np.abs(c));coherent=abs(np.mean(c))
                    ratio=coherent/max(incoherent,1e-20)
                    assert ratio<=1.0000001
                    lines.append(dict(condition=condition,block=block,method=method,harmonic=h,
                               theory_hz=h*fp,time_fft_peak_hz=f[ix],time_fft_peak_amplitude=average[ix],
                               exact_order_coherent_amplitude=coherent,
                               mean_individual_amplitude=incoherent,phasor_retention=ratio))
                    for j in range(3):
                        fusion.append(dict(condition=condition,block=block,method=method,harmonic=h,channel=j+1,
                                           coefficient_real=c[j].real,coefficient_imag=c[j].imag,amplitude=abs(c[j])))
                del y,env,fft,centered
    save_csv(out/'frequency_map.csv',frequency);save_csv(out/'time_envelope_spectra.csv',curves)
    save_csv(out/'fault_harmonics.csv',lines);save_csv(out/'channel_coefficients.csv',fusion)
    summary=[]
    for condition in ['BL','PF50']:
        for method in METHODS:
            for harmonic in range(1,9):
                rr=[r for r in lines if r['condition']==condition and r['method']==method and r['harmonic']==harmonic]
                raw=[r for r in lines if r['condition']==condition and r['method']=='none' and r['harmonic']==harmonic]
                gains=[20*np.log10(max(r['exact_order_coherent_amplitude'],1e-20)/max(b['exact_order_coherent_amplitude'],1e-20)) for r,b in zip(rr,raw)]
                summary.append(dict(condition=condition,method=method,harmonic=harmonic,
                    theory_hz=np.mean([r['theory_hz'] for r in rr]),
                    mean_time_peak_amplitude=np.mean([r['time_fft_peak_amplitude'] for r in rr]),
                    mean_exact_order_amplitude=np.mean([r['exact_order_coherent_amplitude'] for r in rr]),
                    mean_gain_vs_raw_db=np.mean(gains),minimum_gain_db=min(gains),maximum_gain_db=max(gains),
                    mean_phasor_retention=np.mean([r['phasor_retention'] for r in rr])))
    save_csv(out/'harmonic_summary.csv',summary)
    save_json(out/'audit_manifest.json',dict(complete=True,mode='post_hoc_frequency_diagnostic',algorithm_refit=False,
              source_run=str(root),band_hz=[4000,10000],working_fs_hz=fs,window='Hann, coherent-gain corrected',
              frequency_axes='Actual uniformly sampled time-envelope FFT; not a relabelled order axis',
              amplitude='Magnitude of the mean complex channel envelope coefficients; acquisition units',
              exact_order='DFT at h*84/31 on shared carrier-angle grid; same operation for every method',
              source_cfg_geometry_status=cfg['gear_parameter_status'],
              code_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'src_py'/'event_alignment.py']},
              limitations=['One historically inspected record per state; three within-record blocks.',
                  'Spectra and phasor retention cannot alone identify individual fault contacts.',
                  'No SVD or amplitude inverse compensation was added.',
                  'Exact-order gains and time-FFT peak amplitudes are different estimators and are reported separately.']))


def plot_and_report(root,out):
    freq=read_csv(out/'frequency_map.csv');curves=read_csv(out/'time_envelope_spectra.csv')
    lines=read_csv(out/'fault_harmonics.csv');summary=read_csv(out/'harmonic_summary.csv')
    metrics=read_csv(root/'metrics.csv');qa=[]
    fp=pick(freq,'planet_same_side_hz',condition='PF50',block=0)[0]
    fc=pick(freq,'carrier_hz',condition='PF50',block=0)[0]
    fig,axes=plt.subplots(3,1,figsize=(183/25.4,144/25.4))
    fig.subplots_adjust(left=.10,right=.985,bottom=.14,top=.92,hspace=.48)
    maximum=max(np.max(pick(curves,'amplitude',condition='PF50',block=0,method=m)) for m in METHODS)*1.25
    for i,(ax,method) in enumerate(zip(axes,METHODS)):
        f=pick(curves,'frequency_hz',condition='PF50',block=0,method=method)
        a=pick(curves,'amplitude',condition='PF50',block=0,method=method)
        ax.plot(f,a,color=COLORS[i],lw=.9)
        for h in range(1,9):
            ax.axvline(h*fp,color='.65',lw=.5,ls=':',zorder=0)
            if h in [1,2,4,6,8]:ax.text(h*fp,maximum*.88,f'{h}fp',ha='center',fontsize=6,color='.3')
        ax.axvline(2*fc,color='#747F8D',lw=.65,ls='--')
        ax.annotate('2fc',xy=(2*fc,maximum*.7),xytext=(1.0,maximum*.91),fontsize=6,
                    arrowprops={'arrowstyle':'-','lw':.55},ha='center')
        ax.set(xlim=(0,50),ylim=(0,maximum),xlabel='Frequency (Hz)',ylabel='Envelope amplitude')
        mark(ax,chr(97+i),LABEL[i])
        ax.texts[-1].set_x(-.055)
    fig.text(.10,.025,f'PF50, same continuous 27–37 s interval and amplitude scale. Band: 4–10 kHz; 10 s FFT bin spacing: 0.1 Hz.\nfp = fm/31 = {fp:.3f} Hz; fc = fm/84 = {fc:.3f} Hz. Dashed 2fc guide is carrier-related, not the planet-fault fundamental.',fontsize=6,color='.3')
    save(fig,out,'01_fault_frequency_hz',qa)

    fig,axes=plt.subplots(2,2,figsize=(183/25.4,127/25.4))
    fig.subplots_adjust(left=.11,right=.985,bottom=.15,top=.90,wspace=.35,hspace=.6)
    for h,ax in enumerate(axes.flat,1):
        for i,method in enumerate(METHODS):
            f=pick(curves,'frequency_hz',condition='PF50',block=0,method=method)
            a=pick(curves,'amplitude',condition='PF50',block=0,method=method)
            use=np.abs(f-h*fp)<=.8
            ax.plot(f[use],a[use],color=COLORS[i],ls=['-','--','-'][i],lw=[1.6,1.2,.9][i],label=LABEL[i])
        ax.axvline(h*fp,color='.55',ls=':',lw=.65)
        ax.set(xlim=(h*fp-.8,h*fp+.8),ylim=(0,None),xlabel='Frequency (Hz)',ylabel='Envelope amplitude')
        mark(ax,chr(96+h),f'{h}fp = {h*fp:.2f} Hz')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.55,.99),ncol=3)
    fig.text(.11,.025,'Fixed first four fault harmonics, PF50 27–37 s; all methods share the frequency grid and acquisition-amplitude scale.\nVertical guides: predicted frequencies. Subpanels have different y ranges to make small spectral differences visible.',fontsize=6,color='.3')
    save(fig,out,'02_fault_harmonic_zoom',qa)

    fig,axes=plt.subplots(1,2,figsize=(183/25.4,105/25.4))
    fig.subplots_adjust(left=.105,right=.985,bottom=.23,top=.89,wspace=.4)
    for i,method in enumerate(METHODS):
        gain=pick(summary,'mean_gain_vs_raw_db',condition='PF50',method=method)
        retain=pick(summary,'mean_phasor_retention',condition='PF50',method=method)
        if i:axes[0].plot(range(1,9),gain,color=COLORS[i],marker=['o','s'][i-1],ms=3,label=LABEL[i])
        axes[1].plot(range(1,9),retain,color=COLORS[i],marker=['^','o','s'][i],ms=3,label=LABEL[i])
    axes[0].axhline(0,color='.6',lw=.6)
    axes[0].set(xlabel='Planet fault harmonic h',ylabel='Amplitude gain vs unaligned (dB)',xticks=range(1,9))
    axes[0].legend(loc='best',fontsize=6);mark(axes[0],'a','Exact-order harmonic amplitude')
    axes[1].set(xlabel='Planet fault harmonic h',ylabel='Envelope phasor retention',xticks=range(1,9),ylim=(0,1.04))
    axes[1].legend(loc='lower left',fontsize=6);mark(axes[1],'b','Three-channel spectral addition')
    fig.text(.105,.035,'Mean of three PF50 time-block estimates; exact h × 84/31 orders, not time-FFT peak searches.\nRetention = |mean(Cj)| / mean(|Cj|), bounded by 1; it is not a fault-identity score or a new proposed metric.',fontsize=6,color='.3')
    save(fig,out,'03_harmonic_gain_and_retention',qa)

    text=['# 故障特征频率复核：V2 冻结结果','',
          '此前已计算包络阶次谱。本次补充真实 Hz 包络频谱、故障倍频局部放大与精确阶次复系数；未重新拟合算法或改变时移。','',
          '## 结论','',
          '当前 PF50 记录的目标故障频率附近，未对齐信号已存在谱峰。局部包络互相关与新版联合方法只带来小幅谱线变化，联合方法没有表现出稳定的频域优势。第一至第四倍频的频率族/背景指标分别为未对齐 16.26 dB、包络互相关 16.58 dB、联合方法 16.55 dB；联合方法比包络互相关低约 0.03 dB，不能据此宣称诊断谱显著改善。',
          '在本次检查的八条低阶包络谱线上，未对齐三通道的平均复矢量保留率已约为 0.998–1.000，原本就很少发生通道间相消。这限制了仅靠相位对齐恢复这些谱线幅值的空间，但不等于每个故障脉冲都已对齐，更不等于已验证真实故障接触时刻。',
          '因此，本轮结果支持“小幅改变故障频率特征”，不支持“联合方法优于包络互相关”或“大幅增强故障诊断谱”的论文结论。此判断限于已检查记录和当前处理频带。','',
          '## 频率定义','',
          f'PF50 第一测试块平均：行星架 fc={fc:.5f} Hz，啮合 fm={fc*84:.5f} Hz，故障齿同侧接触复现 fp=fm/31={fp:.5f} Hz，二倍频 2fp={2*fp:.5f} Hz。',
          '单个故障齿可分别参与两侧啮合，因此需同时检查 fp 与 2fp 及其倍频；两次响应不必等幅、等间隔，不能预设只有 2fp，也不能把三行星数直接乘到单个故障行星的复现频率上。齿数沿用已有 21/31/84，几何一致性仍待另行核实。',
          '这是 4–10 kHz 加速度响应解调得到的低频包络谱，不是原始加速度的低频啮合谱。最大峰若接近 2fc，应与 fp 区分。真实时间 FFT 长 10 s，频点间距 0.1 Hz，未用补零冒充更高分辨率。','',
          '## 原有频率族/背景指标（三个时间块均值）','',
          '| 方法 | fp、2fp、3fp、4fp 族 dB | 2fp、4fp、6fp、8fp 族 dB |',
          '|---|---:|---:|']
    for method,label in zip(METHODS,['未对齐','局部包络互相关','新版联合方法']):
        a=np.mean(pick(metrics,'family_db',condition='PF50',method=method))
        b=np.mean(pick(metrics,'twice_family_db',condition='PF50',method=method))
        text.append(f'| {label} | {a:.3f} | {b:.3f} |')
    text += ['',
        '以上沿用已有角域谱指标：在每条目标线 ±0.18 行星架阶次范围内取谱峰，四条峰的均方根与指定背景中位数比值取 20log10。两组背景剔除范围不同，不能把两列差值当作两族真实能量之比，也不应称物理信噪比。','',
        '## 各倍频的幅值变化','',
        '| 倍频 | 约 Hz | 包络互相关相对未处理 dB | 联合方法相对未处理 dB | 未处理复矢量保留率 | 联合复矢量保留率 |',
        '|---|---:|---:|---:|---:|---:|']
    for h in range(1,9):
        rows={m:next(r for r in summary if r['condition']=='PF50' and r['method']==m and int(r['harmonic'])==h) for m in METHODS}
        text.append(f'| {h}fp | {float(rows["none"]["theory_hz"]):.2f} | {float(rows["envelope_xcorr"]["mean_gain_vs_raw_db"]):+.3f} | '
                    f'{float(rows["joint_envelope"]["mean_gain_vs_raw_db"]):+.3f} | {float(rows["none"]["mean_phasor_retention"]):.3f} | '
                    f'{float(rows["joint_envelope"]["mean_phasor_retention"]):.3f} |')
    text += ['',
        '此表通过原有共同角坐标，在精确 h×84/31 阶次直接计算每通道复谱系数，先在通道间复数平均再取幅值。增益是三个时间块 dB 比值的均值。它与 Hz 图中局部搜索所得最大谱峰不同，两种估计量没有混用。',
        '复矢量保留率为 |mean(Cj)|/mean(|Cj|)，接近 1 表示这些低阶包络谱线在三通道融合时原本就少有相消；它不意味着每个故障脉冲的时刻都已对齐，也不能替代故障身份验证。低能量线的相位可能不稳定，应与该线幅值共同解释。','',
        '## 分析边界与输出','',
        '- 同一状态只有一条已查看记录，三个时间块不是三个独立试验；不据细小均值差异宣称统计优势。',
        '- 局部包络脉冲变窄，不自动意味着低阶故障谱线显著增强；频域收益应由这里的整段频谱直接检查。',
        '- 未引入 SVD、幅值逆补偿或新的融合权重。本次仍为原时间轴局部时移后的三通道包络均值。',
        '- `frequency_map.csv` 为每块频率；`time_envelope_spectra.csv` 为真实 Hz 曲线；`fault_harmonics.csv` 和 `harmonic_summary.csv` 为逐块/汇总倍频值。',
        '- `channel_coefficients.csv` 为逐通道精确阶次复系数，`audit_manifest.json` 记录算法冻结状态和估计器定义。','']
    (out/'故障频率复核报告.md').write_text('\n'.join(text),encoding='utf-8')
    save_json(out/'figure_qa.json',dict(backend='Python',figures=qa,visual_inspection='pending'))


def main():
    p=argparse.ArgumentParser();p.add_argument('result',type=Path);p.add_argument('--plot-only',action='store_true')
    args=p.parse_args();root=args.result.resolve();out=root/'frequency_audit';out.mkdir(exist_ok=True)
    if not args.plot_only:compute(root,out)
    plot_and_report(root,out)
    print(f'OUTPUT={out}',flush=True)


if __name__=='__main__':main()
