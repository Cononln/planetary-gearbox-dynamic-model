# -*- coding: utf-8 -*-
"""我的 TVMS 与论文图3-2数字化曲线的逐点对比 + 按事件窗口的裂纹统计。"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from pengyue_chapter3_reproduction.parameters import params as P
from pengyue_chapter3_reproduction.tvms.tvms_assembly import TVMSModel

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'tvms_excitation')
ROOT = os.path.join(os.path.dirname(__file__), '..', '..')

eps = P.contact_ratio('sp')
eta = np.linspace(0, 1, 1001)
# 复用已导出的查找表, 避免重积分
lk = np.loadtxt(os.path.join(OUT, 'pair_lookup_sp.csv'), delimiter=',', skiprows=1)
mdl = TVMSModel.__new__(TVMSModel)
mdl.eta = lk[:, 0]; mdl.k_healthy = lk[:, 1]; mdl.k_fault = lk[:, 2]
mdl.epsilon = eps; mdl.fault_tooth = P.fault_tooth

# 裂纹齿啮合窗口: u ∈ [14, 14+ε], θ = (u-14)·2π/zs
u_off = 14.0
theta = np.linspace(0, eps * 2 * np.pi / P.zs, 300)
u = u_off + theta / (2 * np.pi / P.zs)
k_h = mdl.mesh_stiffness(u, healthy_only=True)
k_c = mdl.mesh_stiffness(u)

# 论文数字化曲线
pb = np.load(os.path.join(ROOT, 'fig32_digit_blue.npy'))
pr = np.load(os.path.join(ROOT, 'fig32_digit_red.npy'))
rows = ['theta_rad,k_healthy_mine_N_per_m,k_cracked_mine_N_per_m,'
        'k_healthy_paper_e7,k_cracked_paper_e7']
for i in range(theta.size):
    kh_p = float(np.interp(theta[i], pb[0], pb[1])) * 1e7
    kc_p = float(np.interp(theta[i], pr[0], pr[1])) * 1e7
    rows.append(f'{theta[i]:.6f},{k_h[i]:.6e},{k_c[i]:.6e},{kh_p:.6e},{kc_p:.6e}')
with open(os.path.join(OUT, 'fig3_2_comparison.csv'), 'w') as fh:
    fh.write('\n'.join(rows) + '\n')

# 区间统计 (与论文数字化同一窗口)
def zmean(th, kk, a, b):
    m = (th >= a) & (th <= b)
    return kk[m].mean()

zones = [('single', 0.24, 0.36), ('double_approach', 0.02, 0.20),
         ('double_recess', 0.42, 0.56)]
lines = ['', 'fig 3-2 zone comparison (mine vs paper digitized, N/m):']
for name, a, b in zones:
    mh = zmean(theta, k_h, a, b); mc = zmean(theta, k_c, a, b)
    ph = zmean(pb[0], pb[1] * 1e7, a, b); pc = zmean(pr[0], pr[1] * 1e7, a, b)
    lines.append(f'  {name:15s}: mine healthy {mh:.4e} (paper {ph:.4e}, {100*(mh/ph-1):+.1f}%)'
                 f'   cracked {mc:.4e} (paper {pc:.4e}, {100*(mc/pc-1):+.1f}%)'
                 f'   crack drop mine {100*(1-mc/mh):.2f}% vs paper {100*(1-pc/ph):.2f}%')

# 按事件窗口的统计 (替代被稀释的全程统计)
lines.append('')
lines.append('per-event crack statistics (u in [14, 14+eps]):')
xi_c = u - u_off                       # 裂纹齿对进度 (啮合周期)
msk_s = (xi_c >= eps - 1) & (xi_c < 1)
msk_d = (xi_c < eps - 1) | (xi_c >= 1)
lines.append(f'  single-in-event mean: healthy {k_h[msk_s].mean():.4e}, cracked {k_c[msk_s].mean():.4e}'
             f'  drop {100*(1-k_c[msk_s].mean()/k_h[msk_s].mean()):.2f}%')
lines.append(f'  double-in-event mean: healthy {k_h[msk_d].mean():.4e}, cracked {k_c[msk_d].mean():.4e}'
             f'  drop {100*(1-k_c[msk_d].mean()/k_h[msk_d].mean()):.2f}%')
lines.append(f'  max instantaneous drop {100*(1-(k_c/k_h).min()):.2f}% at theta='
             f'{theta[np.argmin(k_c/k_h)]:.4f} rad')

with open(os.path.join(OUT, 'tvms_validation_report.txt'), 'a') as fh:
    fh.write('\n'.join(lines) + '\n')
print('\n'.join(lines))
