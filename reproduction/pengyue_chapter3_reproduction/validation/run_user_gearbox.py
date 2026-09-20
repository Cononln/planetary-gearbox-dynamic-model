# -*- coding: utf-8 -*-
"""
用户齿轮箱 (21/31/84, 600rpm, 20N·m) 系统响应。

参数来源: 用户提供的参数表 (2026-09-18)。缺失项采用仓库 pg_parameters.m 的
既有估计并登记(见输出 missing_parameters.txt):
  径向支撑 1.2e8/1.5e8/3.0e8/8.0e7 N/m, 扭转 3e3/2e6/0 N·m/rad (太阳/行星/内齿/行星架)
几何: SP 宽安装 a'=39.750056mm, 工作压力角 22.7855°; PR 标准 39.75mm, 20°。
TVMS: 势能法形状, 幅值按用户给定平均啮合刚度标定 (SP 4.5e8, PR 6.0e8 N/m)。
输出: results/user_gearbox/ + figures/user_*.png
"""
import sys, os, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from pengyue_chapter3_reproduction.parameters import params as P
from pengyue_chapter3_reproduction.tvms.external_mesh import external_pair_compliance
from pengyue_chapter3_reproduction.tvms.internal_mesh import internal_pair_compliance
from pengyue_chapter3_reproduction.dynamics.model18 import Model18, rk4_run

matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

ROOT = os.path.join(os.path.dirname(__file__), '..')
FIG = os.path.join(ROOT, 'figures')
RES = os.path.join(ROOT, 'results', 'user_gearbox')
os.makedirs(RES, exist_ok=True)

# ---------------------------------------------------------------- 用户参数
U = dict(zs=21, zp=31, zr=84, Np=3, m=1.5e-3, alpha0=np.deg2rad(20.0),
         alpha_sp=np.deg2rad(22.7855), a_sp=39.750056e-3, a_pr=39.75e-3,
         B=22e-3, B_ring=42e-3, E=207e9, nu=0.29,
         ms=2.00004, mp=0.206247, mr=3.21735, mc=3.28891,
         Is=8.18268e-4, Ip=7.11954e-5, Ir=1.70128e-2, Ic=2.71421e-3,
         kx_s=1.2e8, kx_p=1.5e8, kx_r=3.0e8, kx_c=8.0e7,
         kth_s=3.0e3, kth_r=2.0e6, kth_c=0.0,
         cs=900.0, cp=700.0, cr=1800.0, cc=800.0,
         cth_s=1.0, cth_r=150.0, cth_c=0.5,
         zeta=0.03, e_av=0.5e-6, mu=0.06,
         rpm=600.0, Tout=20.0,
         k_mean_sp=4.5e8, k_mean_pr=6.0e8)

SAVED = {k: getattr(P, k) for k in
         ('zs', 'zp', 'zr', 'Np', 'm_mod', 'alpha0', 'rs', 'rp', 'rr',
          'rb_s', 'rb_p', 'rb_r', 'ra_s', 'ra_p', 'rf_s', 'rf_p', 'ra_r', 'rf_r',
          'rc_carrier', 'E', 'nu', 'G', 'B_face', 'r_int_sun', 'r_int_planet',
          'r_int_ring', 'ms', 'mp_', 'mr', 'mc', 'Is', 'Ip', 'Ir', 'Ic',
          'kx_s', 'kx_p', 'kx_r', 'kx_c', 'kt_s_lin', 'kt_r_lin', 'kt_c_lin',
          'zeta_support', 'zeta_mesh', 'e_av', 'mu_fric', 'Tout',
          'f_sm', 'f_c', 'f_m', 'f_s_fault', 'f_p_fault', 'T_s', 'phi_p')}

def apply_user_params():
    m, a0 = U['m'], U['alpha0']
    P.zs, P.zp, P.zr, P.Np = U['zs'], U['zp'], U['zr'], U['Np']
    P.m_mod = m; P.alpha0 = a0
    P.rs, P.rp, P.rr = m*U['zs']/2, m*U['zp']/2, m*U['zr']/2
    P.rb_s, P.rb_p, P.rb_r = P.rs*np.cos(a0), P.rp*np.cos(a0), P.rr*np.cos(a0)
    P.ra_s, P.rf_s = P.rs + m, P.rs - 1.25*m
    P.ra_p, P.rf_p = P.rp + m, P.rp - 1.25*m
    P.ra_r, P.rf_r = P.rr - m, P.rr + 1.25*m
    P.rc_carrier = U['a_pr']          # 行星架半径 = PR 中心距 39.75mm
    P.E, P.nu = U['E'], U['nu']
    P.G = U['E']/(2*(1+U['nu']))
    P.B_face = U['B']
    P.r_int_sun, P.r_int_planet, P.r_int_ring = P.rf_s/2, P.rf_p/2, P.rf_r/2
    P.ms, P.mp_, P.mr, P.mc = U['ms'], U['mp'], U['mr'], U['mc']
    P.Is, P.Ip, P.Ir, P.Ic = U['Is'], U['Ip'], U['Ir'], U['Ic']
    P.kx_s, P.kx_p, P.kx_r, P.kx_c = U['kx_s'], U['kx_p'], U['kx_r'], U['kx_c']
    P.kt_s_lin = U['kth_s']/P.rb_s**2
    P.kt_r_lin = U['kth_r']/P.rb_r**2
    P.kt_c_lin = 0.0
    P.zeta_support, P.zeta_mesh = U['zeta'], U['zeta']
    P.e_av, P.mu_fric, P.Tout = U['e_av'], U['mu'], U['Tout']
    P.f_sm = U['rpm']/60.0
    P.f_c = P.f_sm*U['zs']/(U['zs']+U['zr'])
    P.f_m = (P.f_sm - P.f_c)*U['zs']
    P.f_s_fault = P.f_m*U['Np']/U['zs']
    P.f_p_fault = P.f_m/U['zp']
    P.T_s = 1.0/P.f_s_fault
    P.phi_p = 2*np.pi*np.arange(U['Np'])/U['Np']

def restore():
    for k, v in SAVED.items():
        setattr(P, k, v)

# ---------------------------------------------------------------- 覆盖后计算
apply_user_params()
ETA = np.linspace(0, 1, 501)
A_SP = P.zs*P.phi_p/(2*np.pi)
A_RP = -P.zr*P.phi_p/(2*np.pi)
SUN_FAULT_TOOTH = 15
CRACK = dict(q=1.5e-3, gamma=np.deg2rad(45.0))
ALPHA_SP_OP = U['alpha_sp']
A_SP_CD = U['a_sp']

SP_S = None

def build_tables():
    global SP_S
    ch = external_pair_compliance(ETA, n_int=300,
                                  alpha_op=ALPHA_SP_OP, a_op=A_SP_CD)
    cc = external_pair_compliance(ETA, crack_sun=CRACK, n_int=300,
                                  alpha_op=ALPHA_SP_OP, a_op=A_SP_CD)
    prh = internal_pair_compliance(ETA, n_int=300)
    SP_S = U['k_mean_sp']/np.mean(1.0/ch['c_total'])
    s_pr = U['k_mean_pr']/np.mean(1.0/prh['c_total'])
    return (ETA, SP_S/ch['c_total']), (ETA, SP_S/cc['c_total']), (ETA, s_pr/prh['c_total'])

dt_glob = (1.0/P.f_m)/512
t0 = time.time()
SP_H, SP_C, PR_H = build_tables()
print('tables built %.1fs' % (time.time()-t0))

EPS_SP = P.contact_ratio('sp')   # 标准 a 的公式对宽安装不适用 → 手工
pb = np.pi*P.m_mod*np.cos(P.alpha0)
ls = np.sqrt(P.ra_s**2-P.rb_s**2); lp = np.sqrt(P.ra_p**2-P.rb_p**2)
lsum = U['a_sp']*np.sin(ALPHA_SP_OP)
EPS_SP = (ls + lp - lsum)/pb
EPS_PR = P.contact_ratio('pr')

def run(case_type, q=None, g=None):
    if case_type == 'healthy':
        m = Model18(SP_H, PR_H, EPS_SP, EPS_PR, A_SP, A_RP,
                    fault_sp_tooth=None, ksp_cracked=None, friction=False,
                    e_av=U['e_av'], alpha_sp=ALPHA_SP_OP, alpha_pr=P.alpha0)
    else:
        m = Model18(SP_H, PR_H, EPS_SP, EPS_PR, A_SP, A_RP,
                    fault_sp_tooth=SUN_FAULT_TOOTH, ksp_cracked=SP_C[1],
                    friction=False, e_av=U['e_av'],
                    alpha_sp=ALPHA_SP_OP, alpha_pr=P.alpha0)
    q0 = static_eq(m)
    dt = (1.0/P.f_m)/512   # 行星轮扭转模态~12.4kHz, 需小步长保稳
    T_TRAN, T_REC = 2.0, 16*P.zs/P.f_m
    T, Q, _ = rk4_run(m, T_TRAN+T_REC, dt, q0=q0, record=True)
    sel = T >= T_TRAN
    ts = T[sel]-T_TRAN; Qs = Q[sel]
    us = m.u_sp[:, 0]
    x = us[0]*(Qs[:,0]-Qs[:,9]) + us[1]*(Qs[:,1]-Qs[:,10]) + m.rbs*Qs[:,2] - m.rbp*Qs[:,11]
    e = U['e_av']*np.sin(2*np.pi*P.f_m*ts + m.te_ph_sp[0])
    F0 = m.Tin/(3*m.rbs)
    kk = np.array([m.k_sp(ts[j]+T_TRAN, 0) for j in range(ts.size)])
    xdyn = x + e - F0/kk
    return ts, x, xdyn, m

def static_eq(m):
    from pengyue_chapter3_reproduction.dynamics.model18 import DOF
    e_keep = m.e_av; m.e_av = 0.0
    q = np.zeros(DOF)
    for _ in range(8):
        g = m.rhs(0.0, q, np.zeros(DOF))*m.M
        J = np.zeros((DOF, DOF))
        for j in range(DOF):
            d = np.zeros(DOF); d[j] = 1e-8
            J[:, j] = (m.rhs(0.0, q+d, np.zeros(DOF))*m.M - g)/d[j]
        dq = np.linalg.solve(J, g)
        q = q - dq
        if np.max(np.abs(dq)) < 1e-12:
            break
    m.e_av = e_keep
    return q

dt_glob = (1.0/P.f_m)/512

def spec(x):
    x = x - x.mean()
    N = x.size
    X = np.fft.rfft(x)/N
    df = 1.0/(N*dt_glob)
    return 2*np.abs(X), df, np.arange(X.size)*df

# ---------------------------------------------------------------- 运行
import csv
CASES = [('healthy', None, None),
         ('sun_q1.5_g45', 1.5, 45.0),
         ('sun_q0.6_g45', 0.6, 45.0),
         ('sun_q1.2_g45', 1.2, 45.0),
         ('sun_q1.8_g45', 1.8, 45.0),
         ('sun_q2.4_g45', 2.4, 45.0),
         ('sun_q1.2_g15', 1.2, 15.0),
         ('sun_q1.2_g30', 1.2, 30.0),
         ('sun_q1.2_g60', 1.2, 60.0),
         ('sun_q1.2_g75', 1.2, 75.0)]
# 裂纹深度比例 = 彭悦: 0.6/1.2/1.8/2.4mm = 齿全高3.375mm的18%/36%/53%/71%
# (彭悦: 2/4/6/8mm / 11.25mm); 角度扫 q=1.2mm = 36% (彭悦 q=4mm=36%)

def run_crack(qmm, gdeg):
    global SP_C
    crack = dict(q=qmm*1e-3, gamma=np.deg2rad(gdeg))
    cc = external_pair_compliance(ETA, crack_sun=crack, n_int=300,
                                  alpha_op=ALPHA_SP_OP, a_op=A_SP_CD)
    sp_c = (ETA, SP_S/ cc['c_total'])
    m = Model18(SP_H, PR_H, EPS_SP, EPS_PR, A_SP, A_RP,
                fault_sp_tooth=SUN_FAULT_TOOTH, ksp_cracked=sp_c[1],
                friction=False, e_av=U['e_av'],
                alpha_sp=ALPHA_SP_OP, alpha_pr=P.alpha0)
    q0 = static_eq(m)
    dt = dt_glob
    T_TRAN, T_REC = 2.0, 16*P.zs/P.f_m
    T, Q, _ = rk4_run(m, T_TRAN+T_REC, dt, q0=q0, record=True)
    sel = T >= T_TRAN
    ts = T[sel]-T_TRAN; Qs = Q[sel]
    us = m.u_sp[:, 0]
    x = us[0]*(Qs[:,0]-Qs[:,9]) + us[1]*(Qs[:,1]-Qs[:,10]) + m.rbs*Qs[:,2] - m.rbp*Qs[:,11]
    e = U['e_av']*np.sin(2*np.pi*P.f_m*ts + m.te_ph_sp[0])
    F0 = m.Tin/(3*m.rbs)
    kk = np.array([m.k_sp(ts[j]+T_TRAN, 0) for j in range(ts.size)])
    xdyn = x + e - F0/kk
    return ts, x, xdyn

results = {}
t0 = time.time()
for nm, q, g in CASES:
    print('---', nm, flush=True)
    if nm == 'healthy':
        ts, xraw, xd, mh = run('healthy')
    else:
        ts, xraw, xd = run_crack(q, g)
    results[nm] = dict(t=ts, raw=xraw, dyn=xd)
    np.savetxt(os.path.join(RES, f'dte_{nm}.csv'),
               np.column_stack([ts, xd*1e6, xraw*1e6]), delimiter=',',
               header='t_s,xsp1dyn_um,xsp1raw_um', comments='')
print('all runs %.0fs' % (time.time()-t0))

def line(x, f):
    xr = x - x.mean(); N = xr.size
    amp = 2*np.abs(np.fft.rfft(xr))/N
    df = 1.0/(N*dt_glob)
    kk = f/df; k0 = int(kk); frac = kk-k0
    return ((1-frac)*amp[k0] + frac*amp[k0+1])*1e6

fm, fs = P.f_m, P.f_s_fault
xh = results['healthy']['dyn']; xc = results['sun_q1.5_g45']['dyn']
print('')
print('谱线 (动态分量, μm):')
print('%-8s %10s %10s' % ('line', 'healthy', 'q1.5g45'))
for nm, f in [('1fm', fm), ('2fm', 2*fm), ('3fm', 3*fm), ('6fm', 6*fm),
              ('fs', fs), ('2fs', 2*fs), ('3fs', 3*fs),
              ('6fm-fs', 6*fm-fs), ('6fm+fs', 6*fm+fs)]:
    print('%-8s %10.5f %10.5f' % (nm, line(xh, f), line(xc, f)))

# 演化 CSV
with open(os.path.join(RES, 'fault_evolution.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['case','line','freq_Hz','amp_um'])
    for nm, q, g in CASES[1:]:
        x = results[nm]['dyn']
        for n in (1, 2, 3, 4):
            w.writerow([nm, f'nfs({n})', n*fs, line(x, n*fs)])
            w.writerow([nm, f'nfm({n})', n*fm, line(x, n*fm)])

# ---------------------------------------------------------------- 图
fig, axs = plt.subplots(2, 1, figsize=(9, 6))
for ax, (x, ttl) in zip(axs, [(xh, 'a) 健康'), (xc, 'b) 太阳轮裂纹 q=1.5mm γ=45°')]):
    tt = results['healthy']['t']
    msk = tt <= 0.5
    ax.plot(tt[msk], x[msk]*1e6, 'b-', lw=0.5)
    ax.set_xlabel('时间 t (s)'); ax.set_ylabel('DTE 动态分量 (μm)'); ax.set_title(ttl)
fig.suptitle('用户齿轮箱(21/31/84, 600rpm, 20N·m) 系统响应时域图')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'user_dte_time.png'), dpi=150); plt.close(fig)

ampH, dfH, frH = spec(xh)
ampC, _, _ = spec(xc)
fig, axs = plt.subplots(2, 2, figsize=(11, 7))
axs[0,0].plot(frH, ampH*1e6, 'b-', lw=0.6); axs[0,0].set_xlim(0, 2000)
axs[0,0].set_title('a) 健康频谱')
axs[0,1].plot(frH, ampC*1e6, 'b-', lw=0.6); axs[0,1].set_xlim(0, 2000)
axs[0,1].set_title('b) 太阳轮裂纹频谱')
msk = frH <= 160
axs[1,0].plot(frH[msk], ampC[msk]*1e6, 'b-', lw=0.8)
for n in range(1, 7):
    axs[1,0].plot(n*fs, line(xc, n*fs), 'r.', ms=6)
axs[1,0].set_title('c) 低频放大 ($nf_s$=24Hz)')
msk = (frH >= 6*fm-120) & (frH <= 6*fm+120)
axs[1,1].plot(frH[msk], ampC[msk]*1e6, 'b-', lw=0.8)
axs[1,1].axvline(6*fm, color='b', lw=0.8)
for sg in (-1, 1):
    for n in (1, 2, 3):
        axs[1,1].plot(6*fm+sg*n*fs, line(xc, 6*fm+sg*n*fs), 'g.', ms=6)
axs[1,1].set_title('d) $6f_m\\pm nf_s$ (6fm=1008Hz)')
for ax in axs.flat: ax.set_xlabel('频率 f (Hz)'); ax.set_ylabel('幅值 (μm)')
fig.suptitle('用户齿轮箱 系统响应频谱图 (动态分量)')
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'user_spectra.png'), dpi=150); plt.close(fig)


# 深度/角度演化分层频谱图
LAYERS = ['#777777', '#5bc2e7', '#2244bb', '#ee8822', '#cc2222']
def layered(names, labels, fname, title, ylab):
    fig = plt.figure(figsize=(11, 6))
    ax = fig.add_subplot(111, projection='3d')
    for j, nm in enumerate(names):
        x = results[nm]['dyn']
        amp, df, fr = spec(x)
        msk = fr <= 2000
        ax.plot(fr[msk], j*2 + 0*fr[msk], amp[msk]*1e6, color=LAYERS[j], lw=0.5)
    ax.set_yticks([0, 2, 4, 6, 8]); ax.set_yticklabels(labels)
    ax.set_ylabel(ylab); ax.set_xlabel('频率 f (Hz)'); ax.set_zlabel('幅值 (μm)')
    ax.set_title(title)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, fname), dpi=150); plt.close(fig)

layered(['healthy','sun_q0.6_g45','sun_q1.2_g45','sun_q1.8_g45','sun_q2.4_g45'],
        ['0','0.6','1.2','1.8','2.4'], 'user_depth_spec.png',
        '用户齿轮箱 不同太阳轮裂纹深度频谱图', '深度 q (mm)')
layered(['sun_q1.2_g15','sun_q1.2_g30','sun_q1.2_g45','sun_q1.2_g60','sun_q1.2_g75'],
        ['15','30','45','60','75'], 'user_angle_spec.png',
        '用户齿轮箱 不同太阳轮裂纹角度频谱图', '角度 γ (°)')

with open(os.path.join(RES, 'missing_parameters.txt'), 'w', encoding='utf-8') as fh:
    fh.write("""用户齿轮箱缺项登记 (2026-09-18)
=========================================
1. 径向支撑刚度(未提供): 采用仓库 pg_parameters.m 既有估计(标注为标定参数):
   太阳轮 1.2e8, 行星轮轴 1.5e8, 内齿圈 3.0e8, 行星架 8.0e7 N/m —— 需确认
2. 扭转支撑刚度(未提供): 同上 3.0e3 / 2.0e6 / 0(行星架) N·m/rad —— 需确认
   (行星架扭转无接地刚度, 通过行星轮-内齿圈回路约束)
3. 支撑阻尼: 用户表给统一 ζ=0.03, 与仓库绝对阻尼(900/700/1800/800 N·s/m)一致, 直接采用
4. 齿侧间隙: 未给 → 线性接触(无间隙)
5. 啮合误差相位: 未给 → Parker 一致相位
6. 裂纹工况: 基准 q=1.5mm γ=45° 齿号15(任意); 如需其他工况请提供
7. 特征频率(推得): fm=168.00, fs=24.00, fp=5.42, fc=2.00 Hz; Ts=1/24 s
8. 几何核验: SP 宽安装 α'=22.7855° 与中心距 39.750056 自洽; ε_SP=1.1502(单齿区85%),
   ε_PR=1.9343; rc=39.75mm; Zr=84≠Zs+2Zp=83 为非标准装配(以工作中心距为准)
9. TVMS: 势能法形状 × 单一全局标定(均值=用户给定 4.5e8/6.0e8 N/m), 未逐工况调参
""")
restore()
print('\nuser gearbox response done -> results/user_gearbox/, figures/user_*.png')
