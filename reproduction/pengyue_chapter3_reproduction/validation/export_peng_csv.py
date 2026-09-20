# -*- coding: utf-8 -*-
"""导出彭悦复现数据为 CSV, 供 pengyue_figures MATLAB 出图 (与 my_gearbox_figures 同格式)。"""
import sys, os, pickle
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from pengyue_chapter3_reproduction.parameters import params as P
from pengyue_chapter3_reproduction.tvms.external_mesh import external_pair_compliance
from pengyue_chapter3_reproduction.tvms.internal_mesh import internal_pair_compliance

ROOT = os.path.join(os.path.dirname(__file__), '..', '..')
RES = os.path.join(ROOT, 'pengyue_chapter3_reproduction', 'results', 'chapter3_sim')
EXC = os.path.join(ROOT, 'pengyue_chapter3_reproduction', 'results', 'tvms_excitation')
DATA = os.path.join(ROOT, 'pengyue_figures', 'data')
os.makedirs(DATA, exist_ok=True)

EPS_SP, EPS_PR = P.contact_ratio('sp'), P.contact_ratio('pr')
FM, FS, FP = P.f_m, P.f_s_fault, P.f_p_fault
DT = 1.0 / (FM * 64)

with open(os.path.join(RES, 'sim_store.pkl'), 'rb') as fh:
    ST = pickle.load(fh)
LK = np.loadtxt(os.path.join(EXC, 'pair_lookup_sp.csv'), delimiter=',', skiprows=1)
LKR = np.loadtxt(os.path.join(EXC, 'pair_lookup_pr.csv'), delimiter=',', skiprows=1)
A_SP = P.zs * P.phi_p / (2 * np.pi)


def k_sp_curve(u, kH, kC=None, fault_tooth=None, z=P.zs, eps=EPS_SP, eta=LK[:, 0]):
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
    return k + np.where(m2, kb, 0.0)


def k_pr_curve(u, kH, kC=None, fault_tooth=None, off=0.0, eps=EPS_PR, eta=LKR[:, 0]):
    n = np.floor(u); xi = u - n
    nf = np.floor(u + off)
    k = np.interp(np.clip(xi, 0, 1) / eps, eta, kH)
    if kC is not None:
        k = np.where(np.mod(nf, P.zp) + 1 == fault_tooth,
                     np.interp(np.clip(xi, 0, 1) / eps, eta, kC), k)
    m2 = xi < eps - 1
    e2 = (1 + xi) / eps
    kb = np.interp(e2, eta, kH)
    if kC is not None:
        kb = np.where(np.mod(nf - 1, P.zp) + 1 == fault_tooth,
                      np.interp(e2, eta, kC), kb)
    return k + np.where(m2, kb, 0.0)


def save(name, cols, hdr):
    np.savetxt(os.path.join(DATA, name), np.column_stack(cols), delimiter=',',
               header=hdr, comments='')
    print('wrote', name)


# ---- 0) 参数表 -----------------------------------------------------------
save('params.csv', [[FM, FS, FP, EPS_SP, EPS_PR, P.zs, P.zp]],
     hdr='fm_Hz,fs_Hz,fp_Hz,eps_sp,eps_pr,zs,zp')

# ---- 1) 响应: 时域(动态分量, 1s) + 频谱(0-2000Hz) + 行星原始DTE ----------
CASES = (['healthy'] +
         ['sun_q%d_g45' % q for q in (2, 3, 4, 6, 8)] +
         ['sun_q4_g%d' % g for g in (15, 30, 45, 60, 75)] +
         ['pl_q%d_g45' % q for q in (2, 3, 4, 6, 8)] +
         ['pl_q4_g%d' % g for g in (15, 30, 45, 60, 75)])
for nm in CASES:
    d = ST[('e0', nm)]
    t, x = d['t'], d['dte_dyn'][0]
    m = t <= 1.0
    save('time_%s.csv' % nm, [t[m], x[m] * 1e6], hdr='t_s,x_um')
    xc = x - x.mean()
    X = np.fft.rfft(xc) / x.size
    fr = np.arange(X.size) * (1.0 / (x.size * DT))
    n2k = int(2000.0 / (fr[1] - fr[0]))
    save('spec_%s.csv' % nm, [fr[:n2k], 2 * np.abs(X[:n2k]) * 1e6],
         hdr='f_Hz,amp_um')

d15 = ST[('e0', 'pl_q3_g45')]
m15 = d15['t'] <= 0.5
save('dte_raw_pl_q3_g45.csv', [d15['t'][m15], d15['dte'][0][m15] * 1e6],
     hdr='t_s,x_um')

# ---- 2) 图3-2: SP 副刚度 vs 转角 (健康 vs 太阳轮裂纹 LK 表) ---------------
u0 = 14.0 + 17.0
th = np.linspace(0, EPS_SP * 2 * np.pi / P.zs, 800)
u = u0 + th / (2 * np.pi / P.zs)
save('tvms_theta_sp.csv', [th, k_sp_curve(u, LK[:, 1]),
     k_sp_curve(u, LK[:, 1], LK[:, 2], P.fault_tooth)],
     hdr='theta_rad,k_healthy_N_per_m,k_cracked_N_per_m')

# ---- 3) 图3-3: 三路啮合相位刚度 (16 周期) ---------------------------------
tc = np.linspace(0, 16.0, 3200)
for i in range(3):
    uu = tc + A_SP[i]
    save('ksptime_healthy_%d.csv' % (i + 1), [tc,
         k_sp_curve(uu, LK[:, 1])], hdr='tTm,k_N_per_m')
    save('ksptime_crack_%d.csv' % (i + 1), [tc,
         k_sp_curve(uu, LK[:, 1], LK[:, 2], P.fault_tooth)],
     hdr='tTm,k_N_per_m')

# ---- 4) 图3-5/3-6 齿对刚度表 (行星轮裂纹 q3 g45) --------------------------
ETA = np.linspace(0, 1, 501)
crack = dict(q=3e-3, gamma=np.deg2rad(45))
chh = external_pair_compliance(ETA, n_int=300)
ch = external_pair_compliance(ETA, crack_planet=crack, n_int=300)
cih = internal_pair_compliance(ETA, n_int=300)
ci = internal_pair_compliance(ETA, crack_planet=crack, n_int=300)
save('tvms_pl_ext.csv', [ETA, 1 / chh['c_total'], 1 / ch['c_total']],
     hdr='eta,k_healthy,k_cracked')
save('tvms_pl_int.csv', [ETA, 1 / cih['c_total'], 1 / ci['c_total']],
     hdr='eta,k_healthy,k_cracked')

# ---- 5) 图3-7/3-8: 37 周期 SP/PR 刚度 (行星轮裂纹) ------------------------
tp = np.linspace(0, 37.0, 7400)
up = FM * tp * 0 + tp          # u 以啮合周期计数
kC_sp_planet = np.interp(LK[:, 0], ETA, 1 / ch['c_total'])
kC_pr_planet = np.interp(LKR[:, 0], ETA, 1 / ci['c_total'])
save('ksp37_healthy.csv', [tp, k_sp_curve(up, LK[:, 1])],
     hdr='tTm,k_N_per_m')
save('ksp37_crack.csv', [tp,
     k_sp_curve(up, LK[:, 1], kC_sp_planet, 15, z=P.zp)], hdr='tTm,k_N_per_m')
save('kpr37_healthy.csv', [tp, k_pr_curve(up, LKR[:, 1])], hdr='tTm,k_N_per_m')
save('kpr37_crack.csv', [tp,
     k_pr_curve(up, LKR[:, 1], kC_pr_planet, fault_tooth=15, off=P.zp / 2.0)],
     hdr='tTm,k_N_per_m')

print('all exports done ->', os.path.abspath(DATA))
