# -*- coding: utf-8 -*-
"""
用户齿轮箱 方案B 数据生成 (误差激励不进入响应, e_av=0)。
谱结构 = 6fm 主导 + nfs 低频族 + 6fm±nfs 边带 (与彭悦图3-10同形态)。
输出到 my_gearbox_figures/data/。
"""
import sys, os, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from pengyue_chapter3_reproduction.parameters import params as P
from pengyue_chapter3_reproduction.tvms.external_mesh import external_pair_compliance
from pengyue_chapter3_reproduction.tvms.internal_mesh import internal_pair_compliance
from pengyue_chapter3_reproduction.dynamics.model18 import Model18, rk4_run, DOF

ROOT = os.path.join(os.path.dirname(__file__), '..')
DATA = os.path.join(ROOT, '..', 'my_gearbox_figures', 'data')
os.makedirs(DATA, exist_ok=True)

U = dict(zs=21, zp=31, zr=84, Np=3, m=1.5e-3, alpha0=np.deg2rad(20.0),
         alpha_sp=np.deg2rad(22.7855), a_sp=39.750056e-3,
         B=22e-3, E=207e9, nu=0.29,
         ms=2.00004, mp=0.206247, mr=3.21735, mc=3.28891,
         Is=8.18268e-4, Ip=7.11954e-5, Ir=1.70128e-2, Ic=2.71421e-3,
         kx_s=1.2e8, kx_p=1.5e8, kx_r=3.0e8, kx_c=8.0e7,
         kth_s=3.0e3, kth_r=2.0e6, kth_c=0.0,
         zeta=0.03, rpm=600.0, Tout=20.0,
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
    P.m_mod = m
    P.alpha0 = a0
    P.rs, P.rp, P.rr = m*U['zs']/2, m*U['zp']/2, m*U['zr']/2
    P.rb_s, P.rb_p, P.rb_r = P.rs*np.cos(a0), P.rp*np.cos(a0), P.rr*np.cos(a0)
    P.ra_s, P.rf_s = P.rs + m, P.rs - 1.25*m
    P.ra_p, P.rf_p = P.rp + m, P.rp - 1.25*m
    P.ra_r, P.rf_r = P.rr - m, P.rr + 1.25*m
    P.rc_carrier = 39.75e-3
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
    P.e_av = 0.0                      # 方案B: 误差不进入响应
    P.mu_fric = 0.06
    P.Tout = U['Tout']
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


apply_user_params()
ETA = np.linspace(0, 1, 501)
A_SP = P.zs*P.phi_p/(2*np.pi)
A_RP = -P.zr*P.phi_p/(2*np.pi)
ALPHA_SP_OP = U['alpha_sp']
A_SP_CD = U['a_sp']
DT = (1.0/P.f_m)/512
T_TRAN, T_REC = 2.0, 16*P.zs/P.f_m
EPS_SP = (np.sqrt(P.ra_s**2 - P.rb_s**2) + np.sqrt(P.ra_p**2 - P.rb_p**2)
          - U['a_sp']*np.sin(ALPHA_SP_OP))/(np.pi*P.m_mod*np.cos(P.alpha0))
EPS_PR = P.contact_ratio('pr')
SUN_TOOTH, PL_TOOTH = 15, 15
F_SUN = P.f_m/P.zs                  # 8 Hz 同相位轮系单路DTE特征频率
F_PL = P.f_m/P.zp                   # 5.42 Hz

print('characteristic: fm=%.2f fs(fam)=%.2f fp=%.2f fc=%.3f epsSP=%.4f epsPR=%.4f'
      % (P.f_m, F_SUN, F_PL, P.f_c, EPS_SP, EPS_PR), flush=True)


def spec(x):
    x = x - x.mean()
    N = x.size
    X = np.fft.rfft(x)/N
    df = 1.0/(N*DT)
    nH = int(2000.0/df)
    return np.arange(nH)*df, 2*np.abs(X[:nH])


def line_spec(f0):
    x = x - x.mean()
    N = x.size
    amp = 2*np.abs(np.fft.rfft(x))/N
    df = 1.0/(N*DT)
    kk = f0/df
    k0 = int(np.floor(kk))
    frac = kk - k0
    return ((1 - frac)*amp[k0] + frac*amp[min(k0 + 1, len(amp) - 1)])*1e6


print('building TVMS tables ...', flush=True)
ch = external_pair_compliance(ETA, n_int=300, alpha_op=ALPHA_SP_OP, a_op=A_SP_CD)
S_SP = U['k_mean_sp']/np.mean(1.0/ch['c_total'])
prh = internal_pair_compliance(ETA, n_int=300)
S_PR = U['k_mean_pr']/np.mean(1.0/prh['c_total'])
KH_SP = S_SP/ch['c_total']


def crack_tables(q_mm, g_deg, member):
    cr = dict(q=q_mm*1e-3, gamma=np.deg2rad(g_deg))
    if member == 'sun':
        cc = external_pair_compliance(ETA, crack_sun=cr, n_int=300,
                                      alpha_op=ALPHA_SP_OP, a_op=A_SP_CD)
        return (S_SP/cc['c_total'],)
    ce = external_pair_compliance(ETA, crack_planet=cr, n_int=300,
                                  alpha_op=ALPHA_SP_OP, a_op=A_SP_CD)
    ci = internal_pair_compliance(ETA, crack_planet=cr, n_int=300)
    return (S_SP/ce['c_total'], S_PR/ci['c_total'])


def ksp_curve(u, kH, kC, fault_tooth, z, eps):
    n = np.floor(u)
    xi = u - n
    k = np.interp(np.clip(xi, 0, 1)/eps, ETA, kH)
    if kC is not None:
        k = np.where(np.mod(n, z) + 1 == fault_tooth,
                     np.interp(np.clip(xi, 0, 1)/eps, ETA, kC), k)
    m2 = xi < eps - 1
    e2 = (1 + xi)/eps
    kb = np.interp(e2, ETA, kH)
    if kC is not None:
        kb = np.where(np.mod(n - 1, z) + 1 == fault_tooth,
                      np.interp(e2, ETA, kC), kb)
    return k + np.where(m2, kb, 0.0)


def kpr_curve(u, kH, kC=None, fault_tooth=None, off=0.0):
    n = np.floor(u)
    nf = np.floor(u + off)          # offset-aware fault cycle (off may be x.5)
    xi = u - n
    k = np.interp(np.clip(xi, 0, 1)/EPS_PR, ETA, kH)
    if kC is not None:
        k = np.where(np.mod(nf, P.zp) + 1 == fault_tooth,
                     np.interp(np.clip(xi, 0, 1)/EPS_PR, ETA, kC), k)
    m2 = xi < EPS_PR - 1
    e2 = (1 + xi)/EPS_PR
    kb = np.interp(e2, ETA, kH)
    if kC is not None:
        kb = np.where(np.mod(nf - 1, P.zp) + 1 == fault_tooth,
                      np.interp(e2, ETA, kC), kb)
    return k + np.where(m2, kb, 0.0)


def build_model(ksp_c=None, sp_fault_gear='sun', fault_pr=None):
    if ksp_c is not None:
        ft = SUN_TOOTH if sp_fault_gear == 'sun' else PL_TOOTH
    else:
        ft = None
    return Model18((ETA, KH_SP), (ETA, S_PR/prh['c_total']),
                   EPS_SP, EPS_PR, A_SP, A_RP,
                   fault_sp_tooth=ft, ksp_cracked=ksp_c,
                   sp_fault_gear=sp_fault_gear, fault_pr=fault_pr,
                   friction=False, e_av=0.0,
                   alpha_sp=ALPHA_SP_OP, alpha_pr=P.alpha0)


def static_eq(m):
    e_keep = m.e_av
    m.e_av = 0.0
    q = np.zeros(DOF)
    for _ in range(8):
        g = m.rhs(0.0, q, np.zeros(DOF))*m.M
        J = np.zeros((DOF, DOF))
        for j in range(DOF):
            d = np.zeros(DOF)
            d[j] = 1e-8
            J[:, j] = (m.rhs(0.0, q + d, np.zeros(DOF))*m.M - g)/d[j]
        dq = np.linalg.solve(J, g)
        q = q - dq
        if np.max(np.abs(dq)) < 1e-12:
            break
    m.e_av = e_keep
    return q


def simulate(m):
    q0 = static_eq(m)
    T, Q, _ = rk4_run(m, T_TRAN + T_REC, DT, q0=q0, record=True)
    sel = T >= T_TRAN
    ts = T[sel] - T_TRAN
    Qs = Q[sel]
    us = m.u_sp[:, 0]
    x = us[0]*(Qs[:, 0] - Qs[:, 9]) + us[1]*(Qs[:, 1] - Qs[:, 10]) \
        + m.rbs*Qs[:, 2] - m.rbp*Qs[:, 11]
    F0 = m.Tin/(3*m.rbs)
    kk = np.array([m.k_sp(t + T_TRAN, 0) for t in ts])
    xd = x - F0/kk
    xpr = (m.u_pr[0, 0]*(Qs[:, 9] - Qs[:, 3]) + m.u_pr[1, 0]*(Qs[:, 10] - Qs[:, 4])
           + m.rbp*Qs[:, 11] - m.rbr*Qs[:, 5])
    kk_pr = np.array([m.k_pr(t + T_TRAN, 0) for t in ts])
    xdpr = xpr - F0/kk_pr
    return ts, xd, xdpr


CASES = (
    [['healthy', None, None, 'sun']] +
    [['sun_q1.5_g45', 1.5, 45.0, 'sun']] +
    [['sun_q%.1f_g45' % q, q, 45.0, 'sun'] for q in (0.4, 0.8, 1.2, 1.6)] +
    [['sun_q1.2_g%d' % g, 1.2, float(g), 'sun'] for g in (15, 30, 60, 75)] +
    [['pl_q1.5_g45', 1.5, 45.0, 'planet']] +
    [['pl_q%.1f_g45' % q, q, 45.0, 'planet'] for q in (0.4, 0.8, 1.2, 1.6)] +
    [['pl_q1.2_g%d' % g, 1.2, float(g), 'planet'] for g in (15, 30, 60, 75)])

results = {}
t0 = time.time()
for nm, q, g, typ in CASES:
    if os.path.exists(os.path.join(DATA, 'spec_%s.csv' % nm)) and        os.path.exists(os.path.join(DATA, 'time_%s.csv' % nm)):
        print('---', nm, 'skip (files exist)', flush=True)
        continue
    print('---', nm, flush=True)
    if q is None:
        m = build_model(None)
    elif typ == 'sun':
        m = build_model(crack_tables(q, g, 'sun')[0])
    else:
        kc_e, kc_i = crack_tables(q, g, 'planet')
        m = build_model(kc_e, sp_fault_gear='planet',
                        fault_pr=dict(planet=0, tooth=PL_TOOTH, k=kc_i,
                                      offset_cycles=P.zp/2.0))
    ts, xd, xdpr = simulate(m)
    results[nm] = dict(t=ts, dyn=xd, dynpr=xdpr)
    f, a = spec(xd)
    np.savetxt(os.path.join(DATA, 'spec_%s.csv' % nm),
               np.column_stack([f, a*1e6]), delimiter=',',
               header='f_Hz,amp_um', comments='')
    msk = ts <= 2.0
    np.savetxt(os.path.join(DATA, 'time_%s.csv' % nm),
               np.column_stack([ts[msk], xd[msk]*1e6]), delimiter=',',
               header='t_s,x_um', comments='')
print('all %d runs %.0fs' % (len(CASES), time.time() - t0), flush=True)

# ---------------------------------------------------------------- TVMS 导出
KC_SP_SUN = crack_tables(1.5, 45.0, 'sun')[0]
th = np.linspace(0, EPS_SP*2*np.pi/P.zs, 600)
kh_pair = np.interp(np.linspace(0, 1, 600), ETA, KH_SP)
kc_pair = np.interp(np.linspace(0, 1, 600), ETA, KC_SP_SUN)
np.savetxt(os.path.join(DATA, 'tvms_theta_sp.csv'),
           np.column_stack([th, kh_pair, kc_pair]), delimiter=',',
           header='theta_rad,k_healthy_N_per_m,k_cracked_N_per_m', comments='')

tp = np.linspace(0, 31.0, 4000)
ce45 = external_pair_compliance(ETA, crack_planet=dict(q=1.5e-3,
                                gamma=np.deg2rad(45.0)), n_int=300,
                                alpha_op=ALPHA_SP_OP, a_op=A_SP_CD)
np.savetxt(os.path.join(DATA, 'ksp31_healthy.csv'),
           np.column_stack([tp, ksp_curve(tp + A_SP[0], KH_SP, None, None,
                                          P.zs, EPS_SP)]),
           delimiter=',', header='tTm,k_N_per_m', comments='')
np.savetxt(os.path.join(DATA, 'ksp31_crack.csv'),
           np.column_stack([tp, ksp_curve(tp + A_SP[0], KH_SP,
                                          S_SP/ce45['c_total'], PL_TOOTH,
                                          P.zp, EPS_SP)]),
           delimiter=',', header='tTm,k_N_per_m', comments='')
ci45 = internal_pair_compliance(ETA, crack_planet=dict(q=1.5e-3,
                                gamma=np.deg2rad(45.0)), n_int=300)
np.savetxt(os.path.join(DATA, 'kpr31_healthy.csv'),
           np.column_stack([tp, kpr_curve(tp, S_PR/prh['c_total'])]),
           delimiter=',', header='tTm,k_N_per_m', comments='')
np.savetxt(os.path.join(DATA, 'kpr31_crack.csv'),
           np.column_stack([tp, kpr_curve(tp, S_PR/prh['c_total'],
                           S_PR/ci45['c_total'], fault_tooth=PL_TOOTH,
                           off=P.zp/2.0)]),
           delimiter=',', header='tTm,k_N_per_m', comments='')

np.savetxt(os.path.join(DATA, 'tvms_pl_ext.csv'),
           np.column_stack([ETA, S_SP/ch['c_total'], S_SP/ce45['c_total']]),
           delimiter=',', header='eta,k_healthy,k_cracked', comments='')
np.savetxt(os.path.join(DATA, 'tvms_pl_int.csv'),
           np.column_stack([ETA, S_PR/prh['c_total'], S_PR/ci45['c_total']]),
           delimiter=',', header='eta,k_healthy,k_cracked', comments='')

# 图3-02: 三行星轮相位刚度 (健康三路相同 = 同相位轮系; 裂纹三路错开7周期)
tc2 = np.linspace(0, 21.0, 3000)
for i in range(3):
    kh_i = ksp_curve(tc2 + A_SP[i], KH_SP, None, None, P.zs, EPS_SP)
    kc_i = ksp_curve(tc2 + A_SP[i], KH_SP, KC_SP_SUN, SUN_TOOTH, P.zs, EPS_SP)
    np.savetxt(os.path.join(DATA, 'ksptime_healthy_%d.csv' % (i + 1)),
               np.column_stack([tc2, kh_i]), delimiter=',',
               header='tTm,k_N_per_m', comments='')
    np.savetxt(os.path.join(DATA, 'ksptime_crack_%d.csv' % (i + 1)),
               np.column_stack([tc2, kc_i]), delimiter=',',
               header='tTm,k_N_per_m', comments='')

# ---------------------------------------------------------------- 演化 CSV
NL = chr(10)
with open(os.path.join(DATA, 'evolution.csv'), 'w') as fh:
    fh.write('case,line,freq_Hz,amp_um' + NL)
    for nm, q, g, typ in CASES[1:]:
        dd = np.loadtxt(os.path.join(DATA, 'spec_%s.csv' % nm), delimiter=',', skiprows=1)
        famp = dd[:, 1]
        df = dd[1, 0] - dd[0, 0]
        fam = F_SUN if typ == 'sun' else F_PL
        def line_spec(f0):
            kk = f0/df
            k0 = int(np.floor(kk)); frac = kk - k0
            return ((1-frac)*famp[k0] + frac*famp[min(k0+1, len(famp)-1)])
        for k in (1, 2, 3, 4):
            fh.write('%s,nf(%d),%.4f,%.6e' % (nm, k, k*fam, line_spec(k*fam)) + NL)
            fh.write('%s,nfm(%d),%.4f,%.6e' % (nm, k, k*P.f_m, line_spec(k*P.f_m)) + NL)
        for sg in (-1, 1):
            for k in (1, 2, 3):
                f0 = 6*P.f_m + sg*k*fam
                tag = '6fm+%d' % k if sg > 0 else '6fm-%d' % k
                fh.write('%s,%s,%.4f,%.6e' % (nm, tag, f0, line_spec(f0)) + NL)

restore()
print('data export done ->', os.path.abspath(DATA))
