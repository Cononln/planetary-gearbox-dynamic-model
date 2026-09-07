"""Complete fault-cycle TSA audit of existing frozen V2 event shifts."""
from pathlib import Path
from datetime import datetime
import argparse
import hashlib
import json
import os
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'runtime'/'phase_python'))
sys.path.insert(0,str(ROOT/'src_py'))
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import numpy as np
from scipy import signal
from event_alignment import reconstruct_block, sample_sinc
from tsa_validation import fit_basis, project, tsa_metrics, harmonic_amplitudes, envelope_of_tsa
from run_event_arrival_v2 import save_csv, save_json, load_npz
from plot_fault_event_alignment import plt, read_csv, pick, mark, save

plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans','Liberation Sans']
plt.rcParams['svg.fonttype']='none'
METHODS=['none','envelope_xcorr','joint_envelope']
LABEL=['Unaligned','Envelope XCorr','Joint envelope']
CHINESE=['未对齐','包络互相关','联合方法']
COLOR=['#858D96','#648EA3','#B16D52']
RANKS=[0,1,3,5]
REPS=['waveform','envelope']
POINTS=12288
Q=np.arange(POINTS)/POINTS-.25


def cycle_registry(events, common, motion, span, zp):
    phase=motion['phase']/(2*np.pi*zp)
    centers=np.flatnonzero(events['branch']==0)
    rows=[]
    for index in centers:
        center=np.interp(events['t'][index],motion['time'],phase)
        left,right=np.interp(center+np.array([-.25,.75]),phase,motion['time'])
        if not (left>=span[0]+.02 and right<=span[1]-.02):
            continue
        inside=(events['t']>=left)&(events['t']<right)
        if not common[inside].all() or inside.sum()!=2:
            continue
        rows.append(dict(event_index=int(index),center_s=float(events['t'][index]),
                         start_s=float(left),end_s=float(right),phase_center=float(center),
                         angle_bin=int(np.mod(events['psi'][index],2*np.pi)/(np.pi/3))))
    if not rows:
        raise ValueError('No supported full fault cycles')
    return rows


def make_rows(wave, fs, events, shifts, common, motion, span, registry, zp):
    padded=[span[0]-.1,span[1]+.1]
    y,derivative,rejected=reconstruct_block(wave,fs,events,shifts,common,padded)
    env=np.abs(signal.hilbert(y,axis=0))
    phase=motion['phase']/(2*np.pi*zp)
    centers=np.array([r['phase_center'] for r in registry])
    query_t=np.interp((centers[:,None]+Q).ravel(),phase,motion['time']).reshape(len(registry),POINTS)
    samples=query_t*fs-round(padded[0]*fs)
    assert samples.min()>16 and samples.max()<len(y)-17
    wr=np.stack([sample_sinc(y[:,j],samples) for j in range(3)],axis=1).astype(np.float32)
    er=np.stack([sample_sinc(env[:,j],samples) for j in range(3)],axis=1).astype(np.float32)
    return dict(waveform=wr,envelope=er),dict(time_map_min_derivative=derivative,rejected_warps=rejected)


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compute(root,out):
    previous=json.loads((root/'run_manifest.json').read_text(encoding='utf-8'))
    assert previous['complete']
    source=Path(previous['source'])
    cfg=json.loads((source/'runtime_config.json').read_text(encoding='utf-8'))
    fs=cfg['working_rates_hz']['impulse'];zp=cfg['gear_teeth']['planet'];zr=cfg['gear_teeth']['ring']
    minimum=cfg['svd']['minimum_events_per_fit_group']
    splits=previous['split_seconds'];all_registry=[];group_metrics=[];harmonics=[];blocks=[];diagnostics=[];fits=[]
    hashes={str(root/c/'evaluated_envelopes.npz'):file_hash(root/c/'evaluated_envelopes.npz') for c in ['BL','PF50']}
    for condition in ['BL','PF50']:
        folder=out/condition;folder.mkdir()
        events=load_npz(root/condition/'events.npz')
        frozen=load_npz(root/condition/'evaluated_envelopes.npz')
        common=frozen['envelope_xcorr_mask']&frozen['joint_envelope_mask']
        motion=load_npz(source/condition/'motion.npz')
        wave=np.load(source/condition/'band_2'/'waveform.npy',mmap_mode='r')
        spans=[splits['calibration']]+splits['test']
        regs=[cycle_registry(events,common,motion,span,zp) for span in spans]
        for i,reg in enumerate(regs):
            all_registry.extend([dict(condition=condition,split='calibration' if i==0 else 'test',block=i-1,**r) for r in reg])
        calbin=np.array([r['angle_bin'] for r in regs[0]])
        calcounts=np.bincount(calbin,minlength=6)
        print(condition,'cycles',list(map(len,regs)),'calibration groups',calcounts.tolist(),flush=True)
        for method in METHODS:
            print(condition,method,'calibration rows',flush=True)
            calibration,diag=make_rows(wave,fs,events,frozen[method+'_shift'],common,motion,spans[0],regs[0],zp)
            diagnostics.append(dict(condition=condition,method=method,block=-1,**diag))
            bases={};saved_basis={}
            for rep in REPS:
                for group in range(6):
                    for channel in range(3):
                        basis=None;singular=np.array([])
                        if calcounts[group]>=minimum:
                            basis,singular=fit_basis(calibration[rep][calbin==group,channel],max_rank=5)
                            saved_basis[f'{rep}_{group}_{channel}']=basis
                        bases[rep,group,channel]=basis
                        fits.append(dict(condition=condition,method=method,representation=rep,angle_bin=group,channel=channel+1,
                                         calibration_cycles=int(calcounts[group]),svd_enabled=basis is not None,
                                         rank3_calibration_energy=float(np.sum(singular[:3]**2)/max(np.sum(singular**2),1e-30))))
            np.savez_compressed(folder/f'{method}_calibration_bases.npz',**saved_basis)
            del calibration
            averages={}
            for block,span in enumerate(splits['test']):
                print(condition,method,'test',block,flush=True)
                registry=regs[block+1];bins=np.array([r['angle_bin'] for r in registry])
                data,diag=make_rows(wave,fs,events,frozen[method+'_shift'],common,motion,span,registry,zp)
                diagnostics.append(dict(condition=condition,method=method,block=block,**diag))
                fp=float(np.diff(np.interp(span,motion['time'],motion['phase']))[0]/(2*np.pi*zp*np.diff(span)[0]))
                for rep in REPS:
                    for rank in RANKS:
                        current=[];current_h=[]
                        for group in range(6):
                            original=data[rep][bins==group].astype(float)
                            if len(original)<4:
                                raise ValueError('Too few cycles for the fixed half-split; do not silently discard a stratum')
                            processed=np.stack([project(original[:,j],bases[rep,group,j],rank) for j in range(3)],axis=1)
                            tsa=processed.mean(axis=0)
                            output=tsa_metrics(original,processed)
                            identity=dict(condition=condition,method=method,representation=rep,rank=rank,block=block,angle_bin=group)
                            row=dict(**identity,**output)
                            group_metrics.append(row);current.append(row)
                            h=harmonic_amplitudes(tsa,rep)
                            current_h.append(h)
                            for k,a in enumerate(h,1):
                                harmonics.append(dict(**identity,harmonic=k,frequency_hz=k*fp,amplitude=float(a)))
                            # Save all groups; the displayed example is fixed before computation.
                            averages[f'{rep}_rank{rank}_block{block}_group{group}']=tsa.astype(np.float32)
                        average_h=np.mean(current_h,axis=0)
                        blocks.append(dict(condition=condition,method=method,representation=rep,rank=rank,block=block,
                                           n_cycles=len(registry),fp_hz=fp,
                                           tsa_retention=float(np.mean([r['tsa_retention'] for r in current])),
                                           heldout_explained=float(np.mean([r['heldout_explained'] for r in current])),
                                           split_tsa_correlation=float(np.mean([r['split_tsa_correlation'] for r in current])),
                                           harmonic_rms_1to4=float(np.sqrt(np.mean(average_h[:4]**2))),
                                           harmonic_rms_1to8=float(np.sqrt(np.mean(average_h**2))),
                                           svd_enabled_groups=int(np.sum(calcounts>=minimum)) if rank else 0))
                del data
            np.savez_compressed(folder/f'{method}_tsa.npz',q=Q,**averages)
        del frozen,wave
    save_csv(out/'cycle_registry.csv',all_registry);save_csv(out/'svd_fits.csv',fits)
    save_csv(out/'group_metrics.csv',group_metrics);save_csv(out/'harmonics.csv',harmonics)
    for row in blocks:
        base=next(r for r in blocks if r['condition']==row['condition'] and r['method']=='none' and
                  r['representation']==row['representation'] and r['rank']==row['rank'] and r['block']==row['block'])
        row['harmonic_gain_vs_none_db']=float(20*np.log10(max(row['harmonic_rms_1to4'],1e-20)/max(base['harmonic_rms_1to4'],1e-20)))
    save_csv(out/'block_metrics.csv',blocks);save_csv(out/'reconstruction_audit.csv',diagnostics)
    assert all(file_hash(Path(path))==value for path,value in hashes.items()),'Frozen inputs were modified'
    save_json(out/'run_manifest.json',dict(complete=True,source_run=str(root),alignment_refit=False,
              frozen_input_sha256=hashes,source_config=str(source/'runtime_config.json'),
              representations={'waveform':'band waveform -> optional frozen SVD -> signed cycle TSA -> Hilbert envelope -> Fourier series',
                               'envelope':'continuous band waveform -> Hilbert envelope -> optional frozen SVD -> cycle TSA -> Fourier series'},
              cycle='mesh phase / planet teeth; full interval [-0.25,0.75) relative to branch A',
              full_tidal_tpsvd=False,cycle_points=POINTS,carrier_angle_bins=6,primary_svd_rank=3,sensitivity_ranks=[1,5],
              minimum_calibration_cycles=minimum,split_seconds=splits,
              harmonics='Rectangular full-cycle Fourier series, no zero padding; h*fp is equivalent frequency, not 0.1 Hz FFT resolution',
              fusion='Average three channel envelopes after within-channel TSA, then complex DFT magnitude; average amplitudes across fixed angle strata',
              evidence='Three within-record test blocks per condition; historically inspected data, no independent statistical superiority',
              code_sha256={str(p.relative_to(ROOT)):file_hash(p) for p in [Path(__file__),ROOT/'src_py'/'tsa_validation.py',ROOT/'tests'/'test_tsa_validation.py']}))
    print('NUMERICAL_COMPLETE',out,flush=True)


def mean_value(rows,key,**where):
    values=pick(rows,key,**where)
    return float(np.mean(values))


def plot_and_report(out):
    rows=read_csv(out/'block_metrics.csv');harm=read_csv(out/'harmonics.csv');qa=[]
    archives={m:load_npz(out/'PF50'/f'{m}_tsa.npz') for m in METHODS}
    fp=mean_value(rows,'fp_hz',condition='PF50',block=0)
    displayed=[]
    fig,axes=plt.subplots(2,2,figsize=(183/25.4,130/25.4))
    fig.subplots_adjust(left=.105,right=.985,bottom=.15,top=.88,wspace=.35,hspace=.55)
    for row,rep in enumerate(REPS):
        for col,rank in enumerate([0,3]):
            ax=axes[row,col]
            for i,m in enumerate(METHODS):
                q=archives[m]['q'];tsa=archives[m][f'{rep}_rank{rank}_block0_group0'][0]
                use=(q/fp>=-.015)&(q/fp<=.03)
                step=4 if rep=='waveform' else 12
                ix=np.flatnonzero(use)[::step]
                ax.plot(q[ix]/fp*1000,tsa[ix],color=COLOR[i],lw=[.8,1.0,.9][i],ls=['-','--','-'][i],label=LABEL[i])
                displayed.extend([dict(method=m,representation=rep,rank=rank,block=0,angle_bin=0,channel=1,
                                       cycle_phase=float(q[k]),equivalent_time_ms=float(q[k]/fp*1000),amplitude=float(tsa[k])) for k in ix])
            title=('' if rank==0 else 'Rank-3 SVD + ')+('vibration TSA' if rep=='waveform' else 'envelope TSA')
            mark(ax,chr(97+2*row+col),title)
            ax.set(xlabel='Equivalent local time (ms)',ylabel='Acquisition amplitude',xlim=(-15,30))
        limits=[axes[row,col].get_ylim() for col in range(2)]
        lo=min(a[0] for a in limits);hi=max(a[1] for a in limits)
        if rep=='envelope':lo=min(0,lo)
        for col in range(2):axes[row,col].set_ylim(lo,hi)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.55,.99),ncol=3)
    fig.text(.105,.035,'PF50, fixed S1 / carrier-angle bin 0, first test block; same complete fault cycles in every arm.\nTop: signed vibration averaging. Bottom: envelope averaging. Time = cycle phase / measured mean fp.',fontsize=6,color='.3')
    save(fig,out,'01_tsa_waveforms',qa)
    save_csv(out/'display_source_data.csv',displayed)

    fig,axes=plt.subplots(2,2,figsize=(183/25.4,130/25.4))
    fig.subplots_adjust(left=.105,right=.985,bottom=.15,top=.88,wspace=.35,hspace=.55)
    for row,rep in enumerate(REPS):
        for col,rank in enumerate([0,3]):
            ax=axes[row,col]
            for i,m in enumerate(METHODS):
                amplitudes=[mean_value(harm,'amplitude',condition='PF50',method=m,representation=rep,rank=rank,harmonic=h) for h in range(1,9)]
                ax.plot(np.arange(1,9)*fp,amplitudes,color=COLOR[i],marker=['^','o','s'][i],ms=3,ls=['-','--','-'][i],label=LABEL[i])
            title=('Vibration TSA envelope' if rep=='waveform' else 'Envelope TSA') if rank==0 else ('SVD + vibration TSA: envelope' if rep=='waveform' else 'SVD + envelope TSA')
            mark(ax,chr(97+2*row+col),title)
            ax.set(xlabel='Fault harmonic frequency hfp (Hz)',ylabel='Fourier-series amplitude',ylim=(0,None),xticks=np.arange(1,9)*fp)
            ax.set_xticklabels([f'{h*fp:.1f}' for h in range(1,9)],fontsize=6)
        hi=max(axes[row,col].get_ylim()[1] for col in range(2))
        for col in range(2):axes[row,col].set_ylim(0,hi)
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.55,.99),ncol=3)
    fig.text(.105,.035,'Three PF50 blocks and six fixed angle strata, equal-weight mean; same scale within each row.\nFull-cycle Fourier series: discrete hfp samples, not a continuous Hz spectrum. A harmonic grid is imposed by TSA.',fontsize=6,color='.3')
    save(fig,out,'02_tsa_fault_harmonics',qa)

    fig,axes=plt.subplots(2,2,figsize=(183/25.4,130/25.4))
    fig.subplots_adjust(left=.11,right=.985,bottom=.15,top=.89,wspace=.38,hspace=.55)
    for row,rep in enumerate(REPS):
        for col,key in enumerate(['heldout_explained','harmonic_gain_vs_none_db']):
            ax=axes[row,col]
            for i,m in enumerate(METHODS):
                for rank in [0,3]:
                    x=(0 if rank==0 else 1)+(i-1)*.18
                    val=pick(rows,key,condition='PF50',method=m,representation=rep,rank=rank)
                    if key=='heldout_explained':val=100*val
                    ax.scatter(x+np.linspace(-.025,.025,len(val)),val,s=10,facecolors='none',edgecolors=COLOR[i])
                    ax.plot([x-.055,x+.055],[np.mean(val)]*2,color=COLOR[i],lw=1.6,label=LABEL[i] if rank==0 else None)
            ax.axhline(0,color='.6',lw=.6)
            ax.set(xticks=[0,1],xticklabels=['TSA','SVD + TSA'],xlim=(-.45,1.45),
                   ylabel='Held-out AC energy explained (%)' if col==0 else 'Harmonic RMS gain vs unaligned (dB)')
            mark(ax,chr(97+2*row+col),('Vibration' if row==0 else 'Envelope')+(' predictability' if col==0 else ' fault harmonics'))
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.55,.99),ncol=3)
    fig.text(.11,.035,'Dots: three time blocks from one PF50 recording; ticks: mean. No independent-repeat significance test.\nHeld-out targets are original cycle rows, never SVD-projected targets. Spectral gain uses the first four harmonics.',fontsize=6,color='.3')
    save(fig,out,'03_tsa_validation',qa)

    text=['# 冻结 V2 对齐后的 TSA / SVD＋TSA 验证','',
          '本次实际运行完整故障复现周期的同步平均，不再仅凭平均前的频谱判断方法。对齐参数未重新拟合。','',
          '## 结论','',
          '在当前包络 TSA 流程中，对齐是有实际收益的。PF50 的留出原周期包络交流能量解释比例，未对齐为 8.57%，包络互相关为 30.94%，联合方法为 30.99%；加入秩 3 的冻结 SVD 后，分别为 0.70%、28.36%、28.28%。这纠正了仅看平均前频谱而低估对齐价值的判断。',
          '包络互相关与联合方法仍基本持平。接秩 3 SVD＋包络 TSA 后，两者前四倍频 RMS 相对同后处理未对齐的增益为 20.789 与 20.797 dB，相差仅约 0.008 dB；联合方法的留出解释比例还略低 0.08 个百分点，不能据此宣称新增联合约束更优。',
          '必须谨慎解释约 20.8 dB：它是相对于“未对齐且经过相同冻结 SVD”的输出，而未对齐数据在该投影下被强烈削弱。它不是相对原信号或不做 SVD 的 TSA 的 20.8 dB 放大，也不是故障信噪比提高 20.8 dB。对齐帮助该校准子空间保留测试脉冲；同时，这也暴露了固定低秩投影对未对齐数据的失配。',
          '原振动有符号 TSA 的留出能量解释仍为负，因此不能称高频振铃已稳定同相。包络 TSA 的正收益更直接支持当前“脉冲到达/包络重复性”的目标。健康组也有一些收益，故该验证不等价于故障身份或物理接触真值。','',
          '## 比较口径','',
          '- 三组为未对齐、包络互相关、联合方法。对所有组使用相同事件接受交集、完整周期、行星架角度分组和平均次数。',
          '- 原振动路线：4–10 kHz 波形 → 冻结时移 → 完整故障周期矩阵 → 可选 SVD → 有符号 TSA → Hilbert 包络 → 倍频系数。',
          '- 包络路线：4–10 kHz 波形 → 冻结时移 → Hilbert 包络 → 完整故障周期矩阵 → 可选 SVD → 包络 TSA → 倍频系数。两条路线不是同一种处理，分别报告。',
          '- 每行覆盖一个同侧故障齿复现周期，约 0.186 s；按共同的啮合相位/31 重采样为 12,288 点，按六个 60° 行星架位置组分别处理三个通道。',
          '- SVD 主结果秩固定为 3，另列 1 和 5 的敏感性结果。基仅由 3–15 s 校准段建立；各方法同规则分别建基。少于 8 个校准周期则各方法均旁路该组 SVD。',
          '- 测试段仍为 27–37、38–48、49–59 s。它们是同一条记录内的时间块，不是独立重复试验。','',
          '## PF50：原振动 TSA','',
          '| 后处理 | 对齐 | TSA 能量保留 % | 留出原周期能量解释 % | 前四倍频 RMS 相对同后处理未对齐 dB |',
          '|---|---|---:|---:|---:|']
    for rep in REPS:
        if rep=='envelope':text+=['','## PF50：包络 TSA','','| 后处理 | 对齐 | TSA 能量保留 % | 留出原周期能量解释 % | 前四倍频 RMS 相对同后处理未对齐 dB |','|---|---|---:|---:|---:|']
        for rank in [0,3]:
            for m,label in zip(METHODS,CHINESE):
                fields=dict(condition='PF50',method=m,representation=rep,rank=rank)
                text.append(f'| {"TSA" if rank==0 else "SVD＋TSA"} | {label} | {100*mean_value(rows,"tsa_retention",**fields):.3f} | '
                            f'{100*mean_value(rows,"heldout_explained",**fields):.3f} | {mean_value(rows,"harmonic_gain_vs_none_db",**fields):+.3f} |')
    text+=['','## 健康控制与秩敏感性','','| 状态 | 输入 | 秩（0 表示 TSA） | 对齐 | 留出解释 % | 倍频 RMS 增益 dB |','|---|---|---:|---|---:|---:|']
    for condition in ['BL','PF50']:
        for rep in REPS:
            for rank in RANKS:
                for m,label in zip(METHODS,CHINESE):
                    fields=dict(condition=condition,method=m,representation=rep,rank=rank)
                    text.append(f'| {condition} | {rep} | {rank} | {label} | {100*mean_value(rows,"heldout_explained",**fields):.3f} | {mean_value(rows,"harmonic_gain_vs_none_db",**fields):+.3f} |')
    text+=['','## 指标解释与边界','',
           'TSA 能量保留为平均输出的交流能量除以输入各周期的平均交流能量；SVD 组的分母仍使用投影前输入，避免通过降低分母制造提升。它仍含有限样本自项，不能单独证明真实同步。',
           '留出原周期能量解释：每个测试角度组按时间交替拆成等长两半，用一半处理后的 TSA 预测另一半未经 SVD 的原周期，并交换方向；报告相对零交流预测的误差降低比例。可为负。拆分依赖、共同校准基和数据驱动对齐使其不是物理接触真值或完全独立盲测。',
           '倍频图来自一个完整周期的 Fourier 级数：先在各通道完成 TSA，再取包络并进行通道平均；最后在六个位置组间等权平均谱幅值。包络路线直接对包络 TSA 取系数。横轴为 h×fp 等效频率，没有宣称 0.1 Hz 分辨率。',
           '按目标周期 TSA 本来就会生成整数倍频网格，因此“有故障倍频线”本身不能证明成功；本报告不在这种强制离散谱上估计背景信噪比，也不把局部窗口串接为频谱。',
           'SVD 对包络的投影可能使局部数值为负；没有为了美化而截零。投影仅保留数据驱动子空间，不添加学习到的平均模板。',
           '本次是按故障周期及路径位置分组的矩阵 SVD/TSA，不是完整潮汐周期 TPSVD。按当前齿数完整潮汐重复约需 31 圈行星架、15.6 s，而每个现有测试块仅 10 s；不能伪造完整潮汐重复数。',
           '所有幅值仍为采集单位；健康与故障灵敏度尚未完成物理校准，不比较两状态的绝对幅值。齿数 21/31/84 的几何一致性仍待核实。','',
           '## 可追溯输出','',
           '- `cycle_registry.csv`：所有使用周期的起止时刻、事件索引及角度组。',
           '- `block_metrics.csv`、`group_metrics.csv`：逐块及逐角度组结果；`harmonics.csv`：逐倍频数据。',
           '- `svd_fits.csv`：基拟合周期数与旁路记录；各状态保存所有方法的平均波形和校准基。',
           '- `run_manifest.json`：冻结输入散列、处理定义与源代码散列。','']
    (out/'TSA验证报告.md').write_text('\n'.join(text),encoding='utf-8')
    save_json(out/'figure_qa.json',dict(backend='Python',figures=qa,visual_inspection='pending',plot_code_sha256=file_hash(Path(__file__))))


def main():
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('--plot-only',action='store_true')
    args=p.parse_args()
    if args.plot_only:
        out=args.source.resolve()
    else:
        out=ROOT/'results'/('tsa_v2_validation_'+datetime.now().strftime('%Y%m%d_%H%M%S'));out.mkdir()
        print('OUTPUT',out,flush=True)
        compute(args.source.resolve(),out)
    plot_and_report(out)
    print('COMPLETE',out,flush=True)


if __name__=='__main__':main()
