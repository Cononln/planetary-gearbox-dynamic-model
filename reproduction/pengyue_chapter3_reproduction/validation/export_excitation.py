# -*- coding: utf-8 -*-
"""
导出彭悦第3章 18DOF 动力学所需的啮合刚度激励。

输出 (results/tvms_excitation/):
  pair_lookup_sp.csv          太阳轮-行星轮单对齿刚度查找表 (健康 + q=3mm,γ=45° 裂纹)
  pair_lookup_pr.csv          行星轮-内齿圈单对齿刚度查找表 (健康)
  sp_pr_tvms_timeseries.csv   三个行星轮 SP/PR 啮合刚度时序 (健康 + 裂纹)
  tvms_harmonics.csv          健康 k_sp1 谐波表 K1..K10
  crack_delta_k_spectrum.csv  裂纹刚度扰动 Δk 频谱 (fs 系列 + fm±n·fs)
  face_width_sensitivity.csv  齿宽 B 敏感性 (B 锁定依据)
  pengyue_tvms_excitation.mat 全部数据 + 标量 (MATLAB 直接加载)
  tvms_validation_report.txt  验证数字汇总
"""
import sys, os, json
import numpy as np
from scipy.io import savemat

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from pengyue_chapter3_reproduction.parameters import params as P
from pengyue_chapter3_reproduction.tvms.external_mesh import external_pair_compliance
from pengyue_chapter3_reproduction.tvms.internal_mesh import internal_pair_compliance
from pengyue_chapter3_reproduction.tvms.tvms_assembly import TVMSModel, sun_fault_event_times

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'tvms_excitation')
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- 1. 单对齿查找表
eta = np.linspace(0, 1, 1001)
crack = dict(q=P.q_crack, gamma=np.deg2rad(P.gamma_deg))
eps_sp = P.contact_ratio('sp')
eps_pr = P.contact_ratio('pr')

comp_sp_h = external_pair_compliance(eta, n_int=400)
comp_sp_c = external_pair_compliance(eta, crack_sun=crack, n_int=400)
comp_pr_h = internal_pair_compliance(eta, n_int=400)

mdl_sp = TVMSModel(comp_sp_h, comp_sp_c, eps_sp, fault_tooth=P.fault_tooth)
mdl_pr = TVMSModel(comp_pr_h, None, eps_pr, fault_tooth=-1)   # 健康内啮合

# ---------------------------------------------------------------- 2. 时序
# 时长 = 16 个太阳轮相对行星架周期 (精确整数个 fm 周期与 fs 周期 → FFT 谱线干净)
n_srel = 16
pts_per_cycle = 256
T_rel = P.zs / P.f_m                      # 0.0890106 s
T = n_srel * T_rel                        # 1.4242 s
N = n_srel * P.zs * pts_per_cycle         # 69632
t = np.arange(N) / (P.f_m * pts_per_cycle)

a_sp = P.zs * P.phi_p / (2 * np.pi)       # [0, 5.6667, 11.3333] (齿对周期)
a_rp = -P.zr * P.phi_p / (2 * np.pi)      # [0, -30.3333, -60.6667]

ksp_h = np.zeros((N, 3)); ksp_c = np.zeros((N, 3)); kpr_h = np.zeros((N, 3))
for i in range(3):
    u_sp = P.f_m * t + a_sp[i]
    ksp_h[:, i] = mdl_sp.mesh_stiffness(u_sp, healthy_only=True)
    ksp_c[:, i] = mdl_sp.mesh_stiffness(u_sp)
    u_pr = P.f_m * t + a_rp[i]
    kpr_h[:, i] = mdl_pr.mesh_stiffness(u_pr, healthy_only=True)

k_mean_sp = ksp_h.mean()
k_mean_pr = kpr_h.mean()
m_eq = 1.0 / (1.0 / P.ms + 1.0 / P.mp_)   # 式(2-43) 等效质量
c_mesh = 2 * P.zeta_mesh * np.sqrt(k_mean_sp * m_eq)   # 啮合阻尼 (ζ=0.07)
omega_n = np.sqrt(k_mean_sp / m_eq)

# ---------------------------------------------------------------- 3. 谐波表
def line_amp(x, f_target, df):
    """精确谱线幅值 (记录长为周期整数, 矩形窗投影)。"""
    n = int(round(f_target / df))
    X = np.fft.rfft(x) / x.size
    return 2 * abs(X[n])

df = 1.0 / T
rows = ['harmonic,frequency_Hz,amplitude_N_per_m,amplitude_over_K1']
K1 = line_amp(ksp_h[:, 0] - ksp_h[:, 0].mean(), P.f_m, df)
for n in range(1, 11):
    f = n * P.f_m
    A = line_amp(ksp_h[:, 0] - ksp_h[:, 0].mean(), f, df)
    rows.append(f'{n}fm,{f:.6f},{A:.6e},{A / K1:.6f}')
with open(os.path.join(OUT, 'tvms_harmonics.csv'), 'w') as fh:
    fh.write('\n'.join(rows) + '\n')

# ---------------------------------------------------------------- 4. Δk 频谱
dk = ksp_h - ksp_c
dk_mean_removed = dk - dk.mean(axis=0)
dk_sum = dk_mean_removed.sum(axis=1)
targets = ([('nfs' if n == 1 else '%dfs' % n, n * P.f_s_fault) for n in (1, 2, 3, 4)] +
           [(f'{n}fm-fs', n * P.f_m - P.f_s_fault) for n in (1, 2, 3, 6)] +
           [(f'{n}fm+fs', n * P.f_m + P.f_s_fault) for n in (1, 2, 3, 6)] +
           [(f'6fm-{m}fs', 6 * P.f_m - m * P.f_s_fault) for m in (2, 3)] +
           [(f'6fm+{m}fs', 6 * P.f_m + m * P.f_s_fault) for m in (2, 3)])
rows = ['line,frequency_Hz,dk1_amp,dk2_amp,dk3_amp,dk_sum_amp_N_per_m']
for name, f in targets:
    a1 = line_amp(dk_mean_removed[:, 0], f, df)
    a2 = line_amp(dk_mean_removed[:, 1], f, df)
    a3 = line_amp(dk_mean_removed[:, 2], f, df)
    aS = line_amp(dk_sum, f, df)
    rows.append(f'{name},{f:.6f},{a1:.6e},{a2:.6e},{a3:.6e},{aS:.6e}')
with open(os.path.join(OUT, 'crack_delta_k_spectrum.csv'), 'w') as fh:
    fh.write('\n'.join(rows) + '\n')

# ---------------------------------------------------------------- 5. 查找表 CSV
# 注意: 模型查找网格是 mdl.eta (n_lookup 点), 必须用它对齐刚度数组
with open(os.path.join(OUT, 'pair_lookup_sp.csv'), 'w') as fh:
    fh.write('eta_norm,k_pair_healthy_N_per_m,k_pair_cracked_N_per_m\n')
    for i in range(mdl_sp.eta.size):
        fh.write(f'{mdl_sp.eta[i]:.9f},{mdl_sp.k_healthy[i]:.6e},{mdl_sp.k_fault[i]:.6e}\n')
with open(os.path.join(OUT, 'pair_lookup_pr.csv'), 'w') as fh:
    fh.write('eta_norm,k_pair_healthy_N_per_m\n')
    for i in range(mdl_pr.eta.size):
        fh.write(f'{mdl_pr.eta[i]:.9f},{mdl_pr.k_healthy[i]:.6e}\n')

# ---------------------------------------------------------------- 6. 时序 CSV
with open(os.path.join(OUT, 'sp_pr_tvms_timeseries.csv'), 'w') as fh:
    fh.write('t_s,k_sp1_healthy,k_sp2_healthy,k_sp3_healthy,'
             'k_sp1_cracked,k_sp2_cracked,k_sp3_cracked,'
             'k_pr1_healthy,k_pr2_healthy,k_pr3_healthy\n')
    for i in range(0, N, 2):          # 隔行写盘减半体积, MAT 里保留全采样
        fh.write(f'{t[i]:.9e},' + ','.join(f'{v:.6e}' for v in
                 [ksp_h[i,0], ksp_h[i,1], ksp_h[i,2],
                  ksp_c[i,0], ksp_c[i,1], ksp_c[i,2],
                  kpr_h[i,0], kpr_h[i,1], kpr_h[i,2]]) + '\n')

# ---------------------------------------------------------------- 7. B 敏感性
rows = ['B_mm,single_zone_mean,double_zone_mean,pair_min,pair_max']
u_test = np.linspace(0, 3, 3001)
xi_test = u_test % 1
msk_s = (xi_test >= eps_sp - 1) & (xi_test < 1)
msk_d = xi_test < eps_sp - 1
B_keep = P.B_face
for B_mm in (4, 6, 8, 10, 12, 20):
    P.B_face = B_mm * 1e-3
    c_h = external_pair_compliance(eta, n_int=250)
    m = TVMSModel(c_h, None, eps_sp)
    k = m.mesh_stiffness(u_test, healthy_only=True)
    rows.append(f'{B_mm},{k[msk_s].mean():.4e},{k[msk_d].mean():.4e},'
                f'{k.min():.4e},{k.max():.4e}')
P.B_face = B_keep
with open(os.path.join(OUT, 'face_width_sensitivity.csv'), 'w') as fh:
    fh.write('\n'.join(rows) + '\n')

# ---------------------------------------------------------------- 8. MAT
savemat(os.path.join(OUT, 'pengyue_tvms_excitation.mat'), {
    't': t,
    'ksp_healthy': ksp_h, 'ksp_cracked': ksp_c, 'kpr_healthy': kpr_h,
    'eta_lookup': eta,
    'kpair_sp_healthy': mdl_sp.k_healthy, 'kpair_sp_cracked': mdl_sp.k_fault,
    'kpair_pr_healthy': mdl_pr.k_healthy,
    'f_m': P.f_m, 'f_s_fault': P.f_s_fault, 'f_p_fault': P.f_p_fault,
    'f_c': P.f_c, 'f_sm': P.f_sm,
    'epsilon_sp': eps_sp, 'epsilon_pr': eps_pr,
    'a_sp_cycles': a_sp, 'a_rp_cycles': a_rp,
    'fault_tooth': P.fault_tooth, 'fault_tooth_offset_from_one': 1,
    'q_crack_m': P.q_crack, 'gamma_crack_deg': P.gamma_deg,
    'k_mean_sp': k_mean_sp, 'k_mean_pr': k_mean_pr,
    'm_eq': m_eq, 'c_mesh': c_mesh, 'omega_n': omega_n,
    'zeta_mesh': P.zeta_mesh, 'zeta_support': P.zeta_support,
    'e_av': P.e_av, 'B_face': P.B_face, 'E': P.E, 'nu': P.nu,
}, do_compression=True)

# ---------------------------------------------------------------- 9. 验证报告
ev = {i + 1: sun_fault_event_times(i, P.fault_tooth, 3) for i in range(3)}
# 单/双齿区掩码: 按行星轮1的局部啮合相位
xi_ts = np.mod(P.f_m * t + a_sp[0], 1.0)
msk_s_ts = (xi_ts >= eps_sp - 1) & (xi_ts < 1)
msk_d_ts = xi_ts < eps_sp - 1
lines = []
A = lines.append
A('Peng Yue Ch.3 TVMS excitation validation (potential-energy method, paper eqs 2-6..2-24, 3-1..3-3)')
A(f'B = {P.B_face*1e3:.1f} mm (locked; see face_width_sensitivity.csv), E = {P.E/1e9:.0f} GPa, nu = {P.nu}')
A(f'epsilon_SP = {eps_sp:.6f}   epsilon_PR = {eps_pr:.6f}')
A(f'f_m = {P.f_m:.6f} Hz   f_s(fault) = {P.f_s_fault:.6f} Hz   f_c = {P.f_c:.6f} Hz')
A(f'T_s (event interval) = {P.T_s:.5f} s  (paper: 0.0297 s)')
A(f'mean k_sp (healthy) = {k_mean_sp:.4e} N/m   mean k_pr = {k_mean_pr:.4e} N/m')
A(f'm_eq = {m_eq:.6f} kg   omega_n = {omega_n:.4f} rad/s   c_mesh = {c_mesh:.2f} N.s/m (zeta=0.07)')
A('')
A('healthy SP zones:')
A(f'  single-zone mean = {ksp_h[:,0][msk_s_ts].mean():.4e} N/m (paper fig 3-2 ~6.2e7)')
A(f'  double-zone mean = {ksp_h[:,0][msk_d_ts].mean():.4e} N/m (paper fig 3-2 peak ~1.19e8, zone mean ~1.1e8)')
A(f'  min = {ksp_h.min():.4e}   max = {ksp_h.max():.4e}')
A('cracked SP zones (q=3mm, gamma=45deg):')
A(f'  single-zone mean = {ksp_c[:,0][msk_s_ts].mean():.4e} ({100*(1-ksp_c[:,0][msk_s_ts].mean()/ksp_h[:,0][msk_s_ts].mean()):.2f}% drop, paper ~2-3%)')
A(f'  double-zone mean = {ksp_c[:,0][msk_d_ts].mean():.4e} ({100*(1-ksp_c[:,0][msk_d_ts].mean()/ksp_h[:,0][msk_d_ts].mean()):.2f}% drop)')
A(f'  max instantaneous drop = {100*(1-(ksp_c/ksp_h).min()):.2f}% (paper fig 3-2 deepest ~5-7%)')
A('')
A('fault event start times per planet (first 3, s):')
for i in range(3):
    A(f'  planet {i+1}: ' + ', '.join(f'{v:.6f}' for v in ev[i+1])
      + f'   (t/T_m = ' + ', '.join(f'{v*P.f_m:.3f}' for v in ev[i+1]) + ')')
A(f'  event interval check: {P.zs/P.f_m/3:.6f} s = T_s/1 (paper 0.0297 s)')
A('')
A('healthy k_sp1 mesh harmonics (K_n/K_1):')
for ln in rows_h if False else open(os.path.join(OUT,'tvms_harmonics.csv')).read().splitlines()[1:]:
    hz, f, amp, ratio = ln.split(',')
    A(f'  {hz:>4s} @ {float(f):8.2f} Hz  K/K1 = {float(ratio):8.5f}')
A('')
A('crack perturbation dk = k_healthy - k_cracked, dominant lines (sum of 3 planets):')
for ln in open(os.path.join(OUT,'crack_delta_k_spectrum.csv')).read().splitlines()[1:]:
    name, f, a1, a2, a3, aS = ln.split(',')
    A(f'  {name:>8s} @ {float(f):9.3f} Hz  dk_sum = {float(aS):.4e} N/m per ... ')
with open(os.path.join(OUT, 'tvms_validation_report.txt'), 'w') as fh:
    fh.write('\n'.join(lines) + '\n')

print('\n'.join(lines[:30]))
print('...')
print('files written to', os.path.abspath(OUT))
