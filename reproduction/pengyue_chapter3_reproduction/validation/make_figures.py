# -*- coding: utf-8 -*-
"""生成彭悦第3章全部信号图 (figures/fig3_*.png) + 演化趋势 CSV。"""
import sys, os, pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from pengyue_chapter3_reproduction.parameters import params as P
from pengyue_chapter3_reproduction.tvms.tvms_assembly import TVMSModel

matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = os.path.join(os.path.dirname(__file__), '..')
FIG = os.path.join(ROOT, 'figures')
RES = os.path.join(ROOT, 'results', 'chapter3_sim')
EXC = os.path.join(ROOT, 'results', 'tvms_excitation')
EPS_SP, EPS_PR = P.contact_ratio('sp'), P.contact_ratio('pr')
fm, fs, fp = P.f_m, P.f_s_fault, P.f_p_fault

with open(os.path.join(RES, 'sim_store.pkl'), 'rb') as fh:
    ST = pickle.load(fh)

LK = np.loadtxt(os.path.join(EXC, 'pair_lookup_sp.csv'), delimiter=',', skiprows=1)
LKR = np.loadtxt(os.path.join(EXC, 'pair_lookup_pr.csv'), delimiter=',', skiprows=1)
A_SP = P.zs * P.phi_p / (2 * np.pi)
LAYERS = ['#777777', '#5bc2e7', '#2244bb', '#ee8822', '#cc2222']

def k_sp_curve(u, kH, kC=None, fault_tooth=None, z=P.zs, eps=EPS_SP, eta=LK[:,0]):
    n = np.floor(u); xi = u - n
    k = np.interp(np.clip(xi, 0, 1) / eps, eta, kH)
    if kC is not None:
        k = np.where(np.mod(n, z) + 1 == fault_tooth,
                     np.interp(np.clip(xi, 0, 1) / eps, eta, kC), k)
    m2 = xi < eps - 1
    e2 = (1 + xi) / eps
    kb = np.interp(e2, eta, kH)
    if kC is not None:
        kb = np.where(np.mod(n - 1, z) + 1 == fault_tooth,
                      np.interp(e2, eta, kC), kb)
    k = k + np.where(m2, kb, 0.0)
    return k

def k_pr_curve(u, kH, kC=None, fault_tooth=None, off=0.0, eps=EPS_PR, eta=LKR[:,0]):
    n = np.floor(u); xi = u - n
    k = np.interp(np.clip(xi, 0, 1) / eps, eta, kH)
    if kC is not None:
        k = np.where(np.mod(n + off, P.zp) + 1 == fault_tooth,
                     np.interp(np.clip(xi, 0, 1) / eps, eta, kC), k)
    m2 = xi < eps - 1
    e2 = (1 + xi) / eps
    kb = np.interp(e2, eta, kH)
    if kC is not None:
        kb = np.where(np.mod(n - 1 + off, P.zp) + 1 == fault_tooth,
                      np.interp(e2, eta, kC), kb)
    k = k + np.where(m2, kb, 0.0)
    return k

# ================================================================ 图 3-2
u0 = 14.0 + 17.0  # 裂纹齿与P1的一次完整啮合窗口
th = np.linspace(0, EPS_SP * 2 * np.pi / P.zs, 800)
u = u0 + th / (2 * np.pi / P.zs)
kh = k_sp_curve(u, LK[:, 1])
kc = k_sp_curve(u, LK[:, 1], LK[:, 2], P.fault_tooth)
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(th, kh / 1e7, 'b-', lw=1.6, label='健康')
ax.plot(th, kc / 1e7, 'r--', lw=1.6, label='含裂纹')
ax.set_xlabel('转角 θ (rad)'); ax.set_ylabel('刚度 $k_{sp1}$ (×10$^7$ N/m)')
ax.legend(); ax.set_title('图3-2 太阳轮-行星轮齿轮副时变啮合刚度 (复现)')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig3_2_tvms.png'), dpi=150); plt.close(fig)

# ================================================================ 图 3-3
tc = np.linspace(0, 16, 4000)   # 啮合周期数
fig, axs = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
for i, ax in enumerate(axs):
    u = tc + A_SP[i]
    kh = k_sp_curve(u, LK[:, 1])
    kc = k_sp_curve(u, LK[:, 1], LK[:, 2], P.fault_tooth)
    ax.plot(tc, kh / 1e7, 'b-', lw=1.0, label='健康')
    ax.plot(tc, kc / 1e7, 'r--', lw=1.0, label='含裂纹')
    ax.set_ylabel(f'$k_{{sp{i+1}}}$ (×10$^7$ N/m)')
axs[0].legend(); axs[-1].set_xlabel('啮合周期 $t/T_m$')
fig.suptitle('图3-3 不同啮合相位下太阳轮-行星轮齿轮副时变啮合刚度 (复现)')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig3_3_mesh_phase.png'), dpi=150); plt.close(fig)

# ============================================ 图 3-5 / 3-6 (行星轮裂纹内外啮合刚度)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from pengyue_chapter3_reproduction.tvms.external_mesh import external_pair_compliance
from pengyue_chapter3_reproduction.tvms.internal_mesh import internal_pair_compliance
crack = dict(q=3e-3, gamma=np.deg2rad(45))
ch = external_pair_compliance(ETA_ := np.linspace(0, 1, 501), crack_planet=crack, n_int=300)
chh = external_pair_compliance(ETA_, n_int=300)
ci = internal_pair_compliance(ETA_, crack_planet=crack, n_int=300)
cih = internal_pair_compliance(ETA_, n_int=300)
for name, comp_h, comp_c, ylab in [
        ('fig3_5_external', chh, ch, '外啮合刚度 (N/m)'),
        ('fig3_6_internal', cih, ci, '内啮合刚度 (N/m)')]:
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(ETA_, 1 / comp_h['c_total'] / 1e7, 'b-', label='健康')
    ax.plot(ETA_, 1 / comp_c['c_total'] / 1e7, 'r--', label='含裂纹')
    ax.set_xlabel('啮合进度 η'); ax.set_ylabel(ylab + ' (×10$^7$ N/m)'); ax.legend()
    ax.set_title(f'图{name.split("_")[1]} 行星轮裂纹(q=3mm,γ=45°) {"外" if "ext" in name else "内"}啮合刚度 (复现)')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, name + '.png'), dpi=150); plt.close(fig)

# ============================================ 图 3-7 / 3-8 (37 啮合周期)
tp = np.linspace(0, 37 / fm, 9000)
up = fm * tp
ksp7_h = k_sp_curve(up, LK[:, 1], kC=None)
# 行星轮裂纹表(501点)插值到查找表网格(2001点)
kC_sp_planet = np.interp(LK[:, 0], ETA_, 1 / ch['c_total'])
ksp7_c = k_sp_curve(up, LK[:, 1], kC_sp_planet, PLANET_FAULT := 15, z=P.zp)
kpr8_h = k_pr_curve(up, LKR[:, 1])
kC_pr_planet = np.interp(LKR[:, 0], ETA_, 1 / ci['c_total'])
kpr8_c = k_pr_curve(up, LKR[:, 1], kC_pr_planet, fault_tooth=15, off=P.zp / 2)
for name, khh, kcc, ylab, ttl in [
        ('fig3_7_sp37', ksp7_h, ksp7_c, '$k_{sp1}$ (×10$^7$ N/m)', '图3-7 太阳轮-行星轮副时变啮合刚度(行星轮裂纹)'),
        ('fig3_8_pr37', kpr8_h, kpr8_c, '$k_{pr1}$ (×10$^7$ N/m)', '图3-8 行星轮-内齿圈副时变啮合刚度(行星轮裂纹)')]:
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(tp, khh / 1e7, 'b-', lw=0.7, label='健康')
    ax.plot(tp, kcc / 1e7, 'r--', lw=0.7, label='含裂纹')
    ax.set_xlabel('啮合周期 $t/T_m$'); ax.set_ylabel(ylab); ax.legend()
    ax.set_title(ttl + ' (复现)')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, name + '.png'), dpi=150); plt.close(fig)

# ================================================================ 响应图工具
def spec(x):
    x = x - x.mean()                     # 去均值(静态变形不入谱)
    n = x.size
    X = np.fft.rfft(x) / n
    df = 1.0 / (n * (1.0 / (P.f_m * 64)))
    return 2 * np.abs(X), df, np.arange(X.size) * df


def dtft(x, f):
    # 任意频率的精确谱线幅值 (DTFT), 避免 fp 不落 FFT 栅格
    x = x - x.mean()
    n = np.arange(x.size)
    return 2 * np.abs(np.dot(x, np.exp(-2j * np.pi * f * n * (1.0/(P.f_m*64))))) / x.size

DTE0 = ST[('e0', 'healthy')]
DTEC = ST[('e0', 'sun_q3_g45')]
D = lambda key: ST[('e0', key)]['dte_dyn']
def line_amp(key, f):
    x = ST[('e0', key)]['dte_dyn'][0]
    x = x - x.mean()
    N = x.size
    amp = 2*np.abs(np.fft.rfft(x))/N
    dfk = 1.0/(N*(1.0/(P.f_m*64)))
    kk = f/dfk
    k0 = int(np.floor(kk)); frac = kk - k0
    return ((1-frac)*amp[k0] + frac*amp[min(k0+1, len(amp)-1)])


# ================================================================ 图 3-9
fig, axs = plt.subplots(2, 1, figsize=(9, 6))
for ax, key, ttl in [(axs[0], 'healthy', 'a) 健康状态下系统响应时域图'),
                     (axs[1], 'sun_q3_g45', 'b) 太阳轮裂纹故障下系统响应时域图')]:
    d = ST[('e0', key)]
    msk = d['t'] <= 0.5
    ax.plot(d['t'][msk], d['dte_dyn'][0][msk] * 1e6, 'b-', lw=0.5)
    ax.set_xlabel('时间 t (s)'); ax.set_ylabel('动态传递误差 DTE (μm)')
    ax.set_title(ttl)
axs[1].annotate('$z_s\\theta_m$', xy=(0.075, axs[1].get_ylim()[1] * 0.9), fontsize=9)
fig.suptitle('图3-9 系统响应时域图 (复现, DTE=动态分量, 扣准静态啮合变形)')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig3_9_dte_time.png'), dpi=150); plt.close(fig)

# ================================================================ 图 3-10
ampH, dfH, frH = spec(DTE0['dte_dyn'][0])
ampC, _, _ = spec(DTEC['dte_dyn'][0])
fig, axs = plt.subplots(2, 2, figsize=(11, 7))
axs[0, 0].plot(frH, ampH * 1e6, 'b-', lw=0.6)
axs[0, 0].set_xlim(0, 2000); axs[0, 0].set_title('a) 健康状态下系统响应频谱图')
axs[0, 1].plot(frH, ampC * 1e6, 'b-', lw=0.6)
axs[0, 1].set_xlim(0, 2000); axs[0, 1].set_title('b) 太阳轮裂纹故障下系统响应频谱图')
msk = frH <= 500
axs[1, 0].plot(frH[msk], ampC[msk] * 1e6, 'b-', lw=0.8)
for n in range(1, 8):
    if n * fs < 500:
        axs[1, 0].plot(n * fs, line_amp('sun_q3_g45', n * fs) * 1e6, 'r.', ms=6)
axs[1, 0].set_title('c) 低频区域频谱放大图 ($nf_s$)')
msk = (frH >= 1000) & (frH <= 1250)
axs[1, 1].plot(frH[msk], ampC[msk] * 1e6, 'b-', lw=0.8)
for n in range(1, 4):
    for sgn in (-1, 1):
        f = 6 * fm + sgn * n * fs
        axs[1, 1].plot(f, line_amp('sun_q3_g45', f) * 1e6, 'g.', ms=6)
axs[1, 1].axvline(6 * fm, color='b', lw=0.8)
axs[1, 1].set_title('d) 高频区域频谱放大图 ($6f_m\\pm nf_s$)')
for ax in axs.flat:
    ax.set_xlabel('频率 f (Hz)'); ax.set_ylabel('幅值 (μm)')
fig.suptitle('图3-10 系统响应频谱图 (复现)')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig3_10_spectra.png'), dpi=150); plt.close(fig)

# ================================================================ 演化图 3-11~3-14
def layered(kind, names, labels, fname, title, ylab=None, mark=None):
    fig = plt.figure(figsize=(11, 6))
    ax = fig.add_subplot(111, projection='3d')
    for j, (nm, lab) in enumerate(zip(names, labels)):
        d = ST[('e0', nm)]
        if kind == 'time':
            msk = d['t'] <= 0.5
            x = d['t'][msk]; y = d['dte_dyn'][0][msk] * 1e6
            ax.plot(x, j * 2 + 0 * x, y, color=LAYERS[j], lw=0.4)
            ax.set_xlabel('时间 t (s)')
        else:
            amp, df, fr = spec(d['dte_dyn'][0])
            msk = fr <= 2000
            ax.plot(fr[msk], j * 2 + 0 * fr[msk], amp[msk] * 1e6, color=LAYERS[j], lw=0.5)
            ax.set_xlabel('频率 f (Hz)')
    ax.set_yticks([0, 2, 4, 6, 8])
    ax.set_yticklabels(labels)
    ax.set_ylabel(ylab or label_axis(kind, names))
    ax.set_zlabel('幅值 (μm)')
    ax.set_title(title + ' (复现)')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, fname), dpi=150); plt.close(fig)

def label_axis(kind, names):
    if 'sun_q' in names[0] and kind == 'spec':
        return '深度 q (mm)'
    if 'sun_q4_g' in names[0]:
        return '角度 γ (°)'
    if 'pl_q' in names[0] and kind == 'spec':
        return '深度 q (mm)'
    return '角度 γ (°)'

sun_depths = ['healthy', 'sun_q2_g45', 'sun_q4_g45', 'sun_q6_g45', 'sun_q8_g45']
sun_angles = ['sun_q4_g15', 'sun_q4_g30', 'sun_q4_g45', 'sun_q4_g60', 'sun_q4_g75']
layered('time', sun_depths, ['0', '2', '4', '6', '8'], 'fig3_11_depth_time.png',
        '图3-11 不同太阳轮裂纹深度下系统响应时域图', ylab='深度 q (mm)')
layered('spec', sun_depths, ['0', '2', '4', '6', '8'], 'fig3_12_depth_spec.png',
        '图3-12 不同太阳轮裂纹深度下系统响应频谱图', ylab='深度 q (mm)')
layered('time', sun_angles, ['15', '30', '45', '60', '75'], 'fig3_13_angle_time.png',
        '图3-13 不同太阳轮裂纹角度下系统响应时域图', ylab='角度 γ (°)')
layered('spec', sun_angles, ['15', '30', '45', '60', '75'], 'fig3_14_angle_spec.png',
        '图3-14 不同太阳轮裂纹角度下系统响应频谱图', ylab='角度 γ (°)')

# ================================================================ 行星轮裂纹 3-15~3-21
DPL = ST[('e0', 'pl_q3_g45')]
fig, ax = plt.subplots(figsize=(9, 3.6))
msk = DPL['t'] <= 0.5
ax.plot(DPL['t'][msk], DPL['dte'][0][msk] * 1e6, 'b-', lw=0.5)
ax.set_xlabel('时间 t (s)'); ax.set_ylabel('DTE (μm)')
ax.set_title('图3-15 行星轮裂纹故障下系统响应时域图 (复现)')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig3_15_pl_time.png'), dpi=150); plt.close(fig)

ampP, dfP, frP = spec(DPL['dte_dyn'][0])
fig, ax = plt.subplots(figsize=(9, 3.6))
ax.plot(frP, ampP * 1e6, 'b-', lw=0.6)
ax.set_xlim(0, 2000); ax.set_xlabel('频率 f (Hz)'); ax.set_ylabel('幅值 (μm)')
ax.set_title('图3-16 行星轮裂纹故障下系统响应频谱图 (复现)')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig3_16_pl_spec.png'), dpi=150); plt.close(fig)

fig, axs = plt.subplots(2, 1, figsize=(9, 6))
msk = frP <= 120
axs[0].plot(frP[msk], ampP[msk] * 1e6, 'b-', lw=0.8)
for n in range(1, 12):
    if n * fp <= 120:
        axs[0].plot(n * fp, line_amp('pl_q3_g45', n * fp) * 1e6, 'r.', ms=6)
axs[0].set_title('a) 低频区域频谱放大图 ($nf_p$)')
msk = (frP >= 1000) & (frP <= 1300)
axs[1].plot(frP[msk], ampP[msk] * 1e6, 'b-', lw=0.8)
for n in range(1, 4):
    for sgn in (-1, 1):
        f = 6 * fm + sgn * n * fp
        axs[1].plot(f, line_amp('pl_q3_g45', f) * 1e6, 'g.', ms=6)
axs[1].axvline(6 * fm, color='b', lw=0.8)
axs[1].set_title('b) 高频区域频谱放大图 ($6f_m\\pm nf_p$)')
for ax in axs: ax.set_xlabel('频率 f (Hz)'); ax.set_ylabel('幅值 (μm)')
fig.suptitle('图3-17 行星轮裂纹故障下系统响应局部频谱图 (复现)')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig3_17_pl_zoom.png'), dpi=150); plt.close(fig)

pl_depths = ['healthy', 'pl_q2_g45', 'pl_q4_g45', 'pl_q6_g45', 'pl_q8_g45']
pl_angles = ['pl_q4_g15', 'pl_q4_g30', 'pl_q4_g45', 'pl_q4_g60', 'pl_q4_g75']
layered('time', pl_depths, ['0', '2', '4', '6', '8'], 'fig3_18_pl_depth_time.png',
        '图3-18 不同行星轮裂纹深度下系统响应时域图', ylab='深度 q (mm)')
layered('spec', pl_depths, ['0', '2', '4', '6', '8'], 'fig3_19_pl_depth_spec.png',
        '图3-19 不同行星轮裂纹深度下系统响应频谱图', ylab='深度 q (mm)')
layered('time', pl_angles, ['15', '30', '45', '60', '75'], 'fig3_20_pl_angle_time.png',
        '图3-20 不同行星轮裂纹角度下系统响应时域图', ylab='角度 γ (°)')
layered('spec', pl_angles, ['15', '30', '45', '60', '75'], 'fig3_21_pl_angle_spec.png',
        '图3-21 不同行星轮裂纹角度下系统响应频谱图', ylab='角度 γ (°)')

# ================================================================ 演化趋势 CSV

rows = ['case,nf_line,freq_Hz,amp_um']
for nm in sun_depths[1:] + sun_angles:
    for n in (1, 2, 3):
        rows.append(f'{nm},nfs({n}),{n*fs:.4f},{line_amp(nm, n*fs)*1e6:.6e}')
        rows.append(f'{nm},nfm({n}),{n*fm:.4f},{line_amp(nm, n*fm)*1e6:.6e}')
for nm in pl_depths[1:] + pl_angles:
    for n in (1, 2, 3):
        rows.append(f'{nm},nfp({n}),{n*fp:.4f},{line_amp(nm, n*fp)*1e6:.6e}')
for sgn in (-1, 1):
    for n in (1, 2, 3):
        rows.append(f'pl_q3_g45,6fm{"+" if sgn>0 else "-"}{n}fp,{6*fm+sgn*n*fp:.4f},'
                    f'{line_amp("pl_q3_g45", 6*fm+sgn*n*fp)*1e6:.6e}')
with open(os.path.join(RES, 'fault_evolution_lines.csv'), 'w') as fh:
    fh.write('\n'.join(rows) + '\n')

print('figures written to', os.path.abspath(FIG))
