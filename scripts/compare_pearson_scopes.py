"""Between-channel Pearson alignment: genuine tidal repeats vs one whole-span lag."""
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
from event_alignment import sample_sinc
from event_arrival import envelope_response,cross_channel_marker_error
from pearson_scope import pearson_lag,pearson,matrix_tsa
from run_event_arrival_v2 import save_csv,save_json,load_npz
from plot_fault_event_alignment import plt,read_csv,pick,mark,save

plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','DejaVu Sans','Liberation Sans']
plt.rcParams['svg.fonttype']='none'
METHODS=['none','whole_envelope','tidal_envelope','whole_waveform','tidal_waveform']
SHOW=METHODS[:3]
LABEL=['Unaligned','Whole-span Pearson','Per-tidal Pearson']
CHINESE=['未对齐','整段固定时移','逐潮汐周期时移']
COLOR=['#858D96','#648EA3','#B16D52']
POINTS=524288
RANKS=[0,1,2]


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def find_periods(motion,zr,zp,start=3.,stop=59.):
    carrier=motion['phase']/(2*np.pi*zr)
    repeat=zp//int(np.gcd(zr,zp))
    origin=float(np.interp(start,motion['time'],carrier))
    count=int(np.floor((np.interp(stop,motion['time'],carrier)-origin)/repeat))
    phase_edges=origin+np.arange(count+1)*repeat
    edges=np.interp(phase_edges,carrier,motion['time'])
    assert count==3,'This protocol expects exactly three full tidal repeats'
    within=np.arange(POINTS)/POINTS*repeat
    query=np.interp((phase_edges[:-1,None]+within).ravel(),carrier,motion['time']).reshape(count,POINTS)
    assert POINTS/np.diff(edges).max()>22000
    return edges,query,repeat


def spectrum(tsa,rep,duration,cycles_per_tidal):
    env=np.abs(signal.hilbert(tsa,axis=-1)).mean(axis=0) if rep=='waveform' else tsa.mean(axis=0)
    amp=2*np.abs(np.fft.rfft(env-env.mean()))/len(env)
    f=np.arange(len(amp))/duration
    fp=cycles_per_tidal/duration;fc=31/duration
    target=fp*np.arange(1,5)
    peaks=np.array([amp[np.abs(f-t)<=.18*fc].max() for t in target])
    use=(f>=fc)&(f<=30*fc)
    for t in target:use &= np.abs(f-t)>.45*fc
    background=float(np.median(amp[use]))
    score=float(20*np.log10(max(np.sqrt(np.mean(peaks*peaks)),1e-20)/max(background,1e-20)))
    h=amp[cycles_per_tidal*np.arange(1,9)]
    keep=f<=60
    return f[keep],amp[keep],h,dict(family_db=score,background=background,
                                  exact_harmonic_rms=float(np.sqrt(np.mean(h[:4]**2))),
                                  fp_hz=fp,fc_hz=fc,fft_grid_hz=1/duration)


def compute(source,out):
    v2=ROOT/'results'/'event_arrival_v2_20260906_072759'
    cfg=json.loads((source/'runtime_config.json').read_text(encoding='utf-8'))
    fs=cfg['working_rates_hz']['impulse'];zr=cfg['gear_teeth']['ring'];zp=cfg['gear_teeth']['planet']
    maximum=round(.008*fs)
    periods=[];fits=[];correlation_curves=[];effective=[];metrics=[];spectra=[];harmonics=[];checks=[];basis_info=[];source_hashes={}
    for condition in ['BL','PF50']:
        folder=out/condition;folder.mkdir()
        motion=load_npz(source/condition/'motion.npz')
        wave=np.load(source/condition/'band_2'/'waveform.npy',mmap_mode='r')
        env=envelope_response(wave,fs)
        edges,query,repeat=find_periods(motion,zr,zp)
        n=len(edges)-1;duration=float(np.mean(np.diff(edges)));cycles=repeat*zr//zp
        for k in range(n):periods.append(dict(condition=condition,tidal=k,start_s=edges[k],end_s=edges[k+1],
                                            duration_s=edges[k+1]-edges[k],carrier_turns=repeat,fault_cycles=cycles))
        print(condition,'full tidal edges',edges.tolist(),flush=True)
        shifts={m:np.zeros((n,3)) for m in METHODS}
        for rep,data in [('envelope',env),('waveform',wave)]:
            for scope in ['whole','tidal']:
                spans=[(edges[0],edges[-1])] if scope=='whole' else list(zip(edges[:-1],edges[1:]))
                for k,(start,end) in enumerate(spans):
                    a,b=round(start*fs),round(end*fs)
                    for j in [1,2]:
                        print(condition,rep,scope,k,'S'+str(j+1),'lag search',flush=True)
                        result,lags,corr=pearson_lag(data[a:b,0],data[a:b,j],maximum)
                        lag=result['lag_samples']/fs
                        if scope=='whole':shifts[f'{scope}_{rep}'][:,j]=lag
                        else:shifts[f'{scope}_{rep}'][k,j]=lag
                        identity=dict(condition=condition,representation=rep,scope=scope,tidal=-1 if scope=='whole' else k,channel=j+1)
                        fits.append(dict(**identity,lag_ms=lag*1000,pearson_zero=result['zero'],pearson_peak=result['peak'],
                                         boundary=result['boundary'],valid=result['valid'],samples=b-a))
                        correlation_curves.extend([dict(**identity,lag_ms=float(l/fs*1000),pearson=float(c)) for l,c in zip(lags,corr)])
        frozen=load_npz(v2/condition/'evaluated_envelopes.npz')
        source_hashes[str(v2/condition/'evaluated_envelopes.npz')]=digest(v2/condition/'evaluated_envelopes.npz')
        markers=dict(offset=frozen['cross_band_marker_offsets'])
        for method in METHODS:
            print(condition,method,'full-tidal resampling and SVD/TSA',flush=True)
            for k in range(n):
                for j in range(3):effective.append(dict(condition=condition,method=method,tidal=k,channel=j+1,lag_ms=shifts[method][k,j]*1000))
                mask=(frozen['t']>=edges[k])&(frozen['t']<edges[k+1])&frozen['main_mask']
                repeated=np.broadcast_to(shifts[method][k],(len(frozen['t']),3))
                spread,pairs=cross_channel_marker_error(markers,repeated,mask)
                checks.append(dict(condition=condition,method=method,tidal=k,cross_band_spread_ms=spread,marker_pairs=pairs))
            saved={}
            for rep,data in [('waveform',wave),('envelope',env)]:
                matrix=np.empty((n,3,POINTS),np.float32)
                for j in range(3):
                    samples=(query+shifts[method][:,j,None])*fs
                    assert samples.min()>16 and samples.max()<len(data)-17
                    matrix[:,j]=sample_sinc(data[:,j],samples).astype(np.float32)
                for k in range(n):
                    for j in [1,2]:checks.append(dict(condition=condition,method=method,tidal=k,
                                                      cross_band_spread_ms=np.nan,marker_pairs=0,
                                                      representation=rep,channel=j+1,pearson=pearson(matrix[k,0],matrix[k,j])))
                for rank in RANKS:
                    means=[]
                    for j in range(3):
                        mean,singular=matrix_tsa(matrix[:,j],rank)
                        means.append(mean)
                        if rank==1:
                            basis_info.append(dict(condition=condition,method=method,representation=rep,channel=j+1,
                                                   rank1_energy=float(singular[0]**2/max(np.sum(singular**2),1e-30)),
                                                   rank2_energy=float(np.sum(singular[:2]**2)/max(np.sum(singular**2),1e-30))))
                    tsa=np.stack(means)
                    saved[f'{rep}_rank{rank}']=tsa.astype(np.float32)
                    f,a,h,row=spectrum(tsa,rep,duration,cycles)
                    identity=dict(condition=condition,method=method,representation=rep,rank=rank)
                    metrics.append(dict(**identity,**row,tidal_periods=n,
                                        tsa_ac_rms=float(np.sqrt(np.mean((tsa-tsa.mean(axis=1,keepdims=True))**2))),
                                        negative_envelope_samples=int((tsa<0).sum()) if rep=='envelope' else 0))
                    spectra.extend([dict(**identity,frequency_hz=float(x),amplitude=float(y)) for x,y in zip(f,a)])
                    harmonics.extend([dict(**identity,harmonic=i+1,frequency_hz=(i+1)*row['fp_hz'],amplitude=float(val)) for i,val in enumerate(h)])
                del matrix
            np.savez_compressed(folder/f'{method}_tsa.npz',**saved)
        del env,wave,frozen
    # Two evaluation types have different fields: keep a rectangular CSV.
    for row in checks:
        row.setdefault('representation','arrival_proxy');row.setdefault('channel',0);row.setdefault('pearson',np.nan)
    save_csv(out/'periods.csv',periods);save_csv(out/'lag_fits.csv',fits);save_csv(out/'lag_curves.csv',correlation_curves)
    save_csv(out/'effective_lags.csv',effective);save_csv(out/'alignment_checks.csv',checks)
    save_csv(out/'metrics.csv',metrics);save_csv(out/'spectra.csv',spectra);save_csv(out/'harmonics.csv',harmonics);save_csv(out/'svd_energy.csv',basis_info)
    assert all(digest(Path(p))==h for p,h in source_hashes.items())
    save_json(out/'run_manifest.json',dict(complete=True,source=str(source),primary_matching_representation='envelope',
              comparison='between-channel S2/S3 to fixed S1, per full tidal period vs one whole matched-span lag',
              alignment_scope='between channels within each tidal repeat; no additional alignment between tidal rows',
              source_fs_hz=fs,tidal_carrier_turns=repeat,tidal_fault_cycles=cycles,cycle_points=POINTS,
              search_half_width_ms=maximum/fs*1000,svd_ranks=RANKS,svd_fit='in-sample three-row matrix; not calibration-frozen or independent heldout',
              downstream='Separate waveform and envelope TSA. Waveform means are enveloped per channel before envelope fusion.',
              frequency_axis='Equivalent-frequency Fourier series of the full angle-resampled tidal mean; not direct uniform-time FFT',
              evaluation='Descriptive same historical recording; optimized Pearson not independent fault-alignment proof',
              original_v2_unchanged_sha256=source_hashes,
              code_sha256={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),ROOT/'src_py'/'pearson_scope.py',ROOT/'tests'/'test_pearson_scope.py']}))
    print('NUMERICAL_COMPLETE',out,flush=True)


def average(rows,key,**where):return float(np.mean(pick(rows,key,**where)))


def plot_report(out):
    lags=read_csv(out/'effective_lags.csv');metrics=read_csv(out/'metrics.csv');spectra=read_csv(out/'spectra.csv');checks=read_csv(out/'alignment_checks.csv');periods=read_csv(out/'periods.csv');fits=read_csv(out/'lag_fits.csv')
    qa=[]
    fig,axes=plt.subplots(1,2,figsize=(183/25.4,91/25.4))
    fig.subplots_adjust(left=.11,right=.985,top=.83,bottom=.24,wspace=.35)
    for j,ax in zip([2,3],axes):
        for i,m in enumerate(SHOW):
            val=pick(lags,'lag_ms',condition='PF50',method=m,channel=j)
            ax.plot([1,2,3],val,marker=['^','o','s'][i],ms=4,ls=['-','--','-'][i],color=COLOR[i],label=LABEL[i])
        ax.axhline(0,color='.7',lw=.5)
        ax.set(xlabel='Full tidal period',ylabel='Shift relative to S1 (ms)',xticks=[1,2,3])
        mark(ax,chr(95+j),f'S{j} relative to S1')
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.55,.99),ncol=3)
    fig.text(.11,.045,'PF50, envelope Pearson search. Three full tidal periods (~15.6 s each); same analyzed span.\nPositive shift: sample moving channel later to match S1. Each tidal period still uses only one scalar shift.',fontsize=6,color='.3')
    save(fig,out,'01_full_tidal_vs_whole_lags',qa)

    fig,axes=plt.subplots(2,2,figsize=(183/25.4,134/25.4))
    fig.subplots_adjust(left=.105,right=.985,top=.88,bottom=.15,hspace=.55,wspace=.35)
    for row,rep in enumerate(['waveform','envelope']):
        for col,rank in enumerate([0,1]):
            ax=axes[row,col]
            for i,m in enumerate(SHOW):
                x=pick(spectra,'frequency_hz',condition='PF50',method=m,representation=rep,rank=rank)
                y=pick(spectra,'amplitude',condition='PF50',method=m,representation=rep,rank=rank)
                ax.plot(x,y,color=COLOR[i],ls=['-','--','-'][i],lw=[1.2,1,.8][i],label=LABEL[i])
            fp=average(metrics,'fp_hz',condition='PF50')
            for h in [1,2,4,6,8]:ax.axvline(h*fp,color='.75',ls=':',lw=.5,zorder=0)
            title=('Vibration TSA envelope' if rep=='waveform' else 'Envelope TSA') if rank==0 else ('Rank-1 SVD + vibration TSA' if rep=='waveform' else 'Rank-1 SVD + envelope TSA')
            mark(ax,chr(97+2*row+col),title)
            ax.set(xlim=(0,50),ylim=(0,None),xlabel='Equivalent frequency (Hz)',ylabel='Envelope amplitude')
        hi=max(axes[row,c].get_ylim()[1] for c in range(2))
        for c in range(2):axes[row,c].set_ylim(0,hi)
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.55,.99),ncol=3)
    fig.text(.105,.035,'PF50: mean of three full tidal rows; envelope-based channel lag estimation in all compared arms.\nDotted lines: 1, 2, 4, 6, 8 times fp. Shared y scale within each row; Fourier-series grid ~0.064 Hz.',fontsize=6,color='.3')
    save(fig,out,'02_full_tidal_tsa_spectra',qa)

    fig,axes=plt.subplots(1,2,figsize=(183/25.4,98/25.4))
    fig.subplots_adjust(left=.11,right=.985,top=.83,bottom=.26,wspace=.4)
    for i,m in enumerate(SHOW):
        pcc=np.array([average(checks,'pearson',condition='PF50',method=m,tidal=k,representation='envelope') for k in range(3)])
        err=pick(checks,'cross_band_spread_ms',condition='PF50',method=m,representation='arrival_proxy')
        for ax,val in zip(axes,[pcc,err]):
            ax.scatter(i+np.linspace(-.04,.04,3),val,color=COLOR[i],s=12)
            ax.plot([i-.12,i+.12],[np.mean(val)]*2,color=COLOR[i],lw=1.5)
    for ax in axes:ax.set(xticks=[0,1,2],xticklabels=['Unaligned','Whole span','Per tidal'],xlim=(-.45,2.45))
    axes[0].set(ylabel='Within-period envelope Pearson',ylim=(-.05,1.05));mark(axes[0],'a','Optimized-signal similarity')
    axes[1].set(ylabel='Cross-band marker spread (ms)',ylim=(0,None));mark(axes[1],'b','Fault-event response proxy')
    fig.text(.11,.05,'Dots: three dependent full tidal periods, not independent experimental repeats.\nPearson is in-sample; cross-band markers are an arrival proxy, not contact truth. Matching uses the envelope.',fontsize=6,color='.3')
    save(fig,out,'03_alignment_similarity_and_events',qa)

    manifest=json.loads((out/'run_manifest.json').read_text(encoding='utf-8'))
    original=Path(manifest['source']);motion=load_npz(original/'PF50'/'motion.npz')
    events=load_npz(ROOT/'results'/'event_arrival_v2_20260906_072759'/'PF50'/'events.npz')
    start=float(next(r['start_s'] for r in periods if r['condition']=='PF50' and r['tidal']=='0'))
    duration=average(periods,'duration_s',condition='PF50')
    event_index=np.flatnonzero((events['branch']==0)&(events['t']>=start+.2))[0]
    phi=np.interp([start,events['t'][event_index]],motion['time'],motion['phase'])/(2*np.pi*84)
    center=float((phi[1]-phi[0])/31*duration)
    time=np.arange(POINTS)/POINTS*duration
    selected=np.flatnonzero((time-center>=-.02)&(time-center<=.04))
    arrays={m:load_npz(out/'PF50'/f'{m}_tsa.npz') for m in SHOW};display=[]
    fig,axes=plt.subplots(2,2,figsize=(183/25.4,128/25.4))
    fig.subplots_adjust(left=.105,right=.985,top=.88,bottom=.15,hspace=.55,wspace=.35)
    for row,rep in enumerate(['waveform','envelope']):
        for col,rank in enumerate([0,1]):
            ax=axes[row,col]
            for i,m in enumerate(SHOW):
                # Show S2, which is actually shifted; S1 is the unchanged reference.
                yy=arrays[m][f'{rep}_rank{rank}'][1,selected]
                xx=(time[selected]-center)*1000
                ax.plot(xx,yy,color=COLOR[i],ls=['-','--','-'][i],lw=[1.1,1,.8][i])
                display.extend([dict(method=m,representation=rep,rank=rank,channel=2,event_index=int(event_index),
                                     equivalent_local_ms=float(x),amplitude=float(y)) for x,y in zip(xx,yy)])
            title=('Vibration TSA' if rep=='waveform' else 'Envelope TSA') if rank==0 else ('SVD + vibration TSA' if rep=='waveform' else 'SVD + envelope TSA')
            mark(ax,chr(97+2*row+col),title)
            ax.set(xlabel='Equivalent local time (ms)',ylabel='S2 acquisition amplitude',xlim=(-20,40))
        lo=min(axes[row,c].get_ylim()[0] for c in range(2));hi=max(axes[row,c].get_ylim()[1] for c in range(2))
        if rep=='envelope':lo=min(lo,0)
        for c in range(2):axes[row,c].set_ylim(lo,hi)
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.55,.99),ncol=3)
    fig.text(.105,.035,'S2, three-tidal-period TSA; fixed first branch-A candidate after 0.2 s in the first tidal period.\nAll methods use the same window and y scale. The candidate is a reference location, not contact ground truth.',fontsize=6,color='.3')
    save(fig,out,'04_full_tidal_tsa_waveform',qa)
    save_csv(out/'waveform_display_data.csv',display)

    text=['# 潮汐周期内通道对齐与整段通道对齐：Pearson 互相关比较','',
          '本轮按“每个潮汐周期内 S2、S3 分别对齐 S1”实现，与“整段 S2、S3 各用一个固定时移”比较；不是同一通道的不同潮汐周期相互对齐。','',
          '## 结论','',
          '在当前 PF50 记录、±8 ms 的搜索范围和固定 S1 参考下，两种对齐尺度的结果基本一致。整段法的 S2、S3 时移分别约为 +0.03934、−0.00991 ms；逐潮汐周期估计的变化仅为微秒量级。该数值精度来自亚采样插值，不代表已经标定出同等精度的物理传播时延。',
          '包络 Pearson 从未对齐 0.77923 变为整段 0.78072、逐潮汐周期 0.78076；跨频带到达标记离散从 2.02935 ms 变为 2.03355、2.03342 ms，没有观察到故障候选脉冲到达一致性的改善。',
          '接相同秩 1 的完整潮汐矩阵 SVD＋包络 TSA 后，目标频率族/背景比分别为 25.171、25.174、25.177 dB；差异不足以支持逐潮汐周期更优。改用秩 2 也没有稳定收益。振动波形相关的对照结果同样没有跨后处理设置的一致优势。',
          '该结果不否定潮汐周期用于组织完整故障信息的意义；它说明每个约 15.6 s 潮汐周期只给一个通道时移，仍是在较长时间内拟合整体相似性，不能据此声称已跟踪周期内随行星架转动的故障脉冲偏移。','',
          '## 处理定义','',
          '- Pearson＋互相关：逐滞后计算去均值、局部方差归一化的 Pearson 相关系数，取最大正相关，抛物线细化时移。不是先用 Pearson 再额外叠加另一个独立算法。',
          '- 主比较在 4–10 kHz 包络上估计时移，平滑尺度沿用 0.12 ms；另完整保留振动波形上的 Pearson 搜索作为对照。',
          '- 参考传感器固定为 S1；两种范围均搜索约 ±8 ms，不按结果更换参考、不反转极性、不循环回卷。',
          '- 潮汐周期采用当前完整状态复现定义：31 圈行星架、84 个同侧故障复现周期。起点 3 s，按已有振动相位确定边界，仅使用三个完整潮汐周期；整段法使用同一完整范围，双方均不使用末尾残段。',
          '- 每条潮汐行重采样为 524,288 点；这是角坐标网格，原始工作采样率仍为 51,200 Hz。按通道分别做三行矩阵的 SVD 与 TSA，统一比较不做 SVD、秩 1、秩 2。',
          '- 此处 SVD 在三行数据上直接做描述性降噪；不是此前局部事件验证中的校准冻结 SVD，不能将两套结果混合或宣称独立盲测。','',
          '## 完整周期边界','',
          '| 状态 | 周期 | 起点 s | 终点 s | 时长 s |','|---|---:|---:|---:|---:|']
    for r in periods:text.append(f'| {r["condition"]} | {int(r["tidal"])+1} | {float(r["start_s"]):.5f} | {float(r["end_s"]):.5f} | {float(r["duration_s"]):.5f} |')
    text+=['','## PF50 包络 Pearson 估计的时移','','| 方式 | 周期 | S2→S1 ms | S3→S1 ms |','|---|---:|---:|---:|']
    for m,label in zip(SHOW[1:],CHINESE[1:]):
        for k in range(3):
            text.append(f'| {label} | {k+1} | {average(lags,"lag_ms",condition="PF50",method=m,tidal=k,channel=2):+.5f} | {average(lags,"lag_ms",condition="PF50",method=m,tidal=k,channel=3):+.5f} |')
    text+=['','## 包络匹配：相似性与事件时刻代理','','| 状态 | 方式 | 周期内包络 Pearson | 跨频带到达标记离散 ms |','|---|---|---:|---:|']
    for condition in ['BL','PF50']:
        for m,label in zip(SHOW,CHINESE):
            text.append(f'| {condition} | {label} | {average(checks,"pearson",condition=condition,method=m,representation="envelope"):.5f} | {average(checks,"cross_band_spread_ms",condition=condition,method=m,representation="arrival_proxy"):.5f} |')
    text+=['','## 潮汐 TSA 后的目标频率族/背景比','','| 状态 | 时移估计 | 方式 | TSA 输入 | SVD 秩（0=不做） | 频率族/背景 dB | 精确前四倍频 RMS |','|---|---|---|---|---:|---:|---:|']
    for r in metrics:
        m=r['method'];match='无' if m=='none' else m.split('_',1)[1]
        scope='未对齐' if m=='none' else ('整段' if m.startswith('whole') else '逐潮汐周期')
        text.append(f'| {r["condition"]} | {match} | {scope} | {r["representation"]} | {r["rank"]} | {float(r["family_db"]):.3f} | {float(r["exact_harmonic_rms"]):.5f} |')
    boundary=sum(r['boundary']=='True' for r in fits)
    text+=['','## 解释边界','',
           f'- 搜索边界命中 {boundary}/{len(fits)} 次。若命中边界，不能解释为已经找到无约束的最佳时移。',
           '- 一个完整潮汐周期仍约 15.6 s；每周期一个时移，无法跟踪该周期内约 31 圈行星架带来的路径变化。这一试验检验分段尺度，不等价于逐故障事件对齐。',
           '- 周期内 Pearson 是经过优化后的同信号相似性，不能独立证明故障脉冲对齐。跨频带标记也只是测量响应代理，不是真实接触时间。',
           '- 频谱是完整潮汐平均波形的 Fourier 级数，等效频率网格约 0.064 Hz；不是将事件短窗串接，也不是原时间轴直接 FFT。',
           '- 频率族/背景比沿用目标线邻域 ±0.18 行星架频率取峰、背景剔除 ±0.45 行星架频率的规则，不是物理信噪比。',
           '- 三个周期属于同一条已查看记录，只做描述性比较，不作独立重复显著性声明。健康与故障仅比较各自处理收益，不比较未校准的绝对物理幅值。',
           '- 当前齿数 21/31/84 的几何一致性仍待核实；此次保持既有定义，不另行改变齿轮模型。','',
           '## 文件','',
           '`periods.csv` 保存真实完整周期边界；`lag_fits.csv`、`lag_curves.csv` 保存时移与完整搜索曲线；`alignment_checks.csv` 保存逐周期检查；`spectra.csv`、`harmonics.csv` 与 `metrics.csv` 保存全部频域结果。各状态目录保留完整三通道 TSA 数组，图以相同尺度导出。','']
    (out/'潮汐周期与整段对齐比较报告.md').write_text('\n'.join(text),encoding='utf-8')
    save_json(out/'figure_qa.json',dict(backend='Python',figures=qa,visual_inspection='pending',plot_code_sha256=digest(Path(__file__))))


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=ROOT/'results'/'event_alignment_v1_20260906_053745');p.add_argument('--plot-only',type=Path)
    args=p.parse_args()
    if args.plot_only:out=args.plot_only.resolve()
    else:
        out=ROOT/'results'/('pearson_tidal_vs_whole_'+datetime.now().strftime('%Y%m%d_%H%M%S'));out.mkdir()
        print('OUTPUT',out,flush=True);compute(args.source.resolve(),out)
    plot_report(out);print('COMPLETE',out,flush=True)


if __name__=='__main__':main()
