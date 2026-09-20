# -*- coding: utf-8 -*-
"""
Chaari/Sainsot 齿基(轮体)柔度 —— 论文式(2-22), 文献[71]:
  Chaari F, Fakhfakh T, Haddar M. Analytical modelling of spur gear tooth crack
  and influence on gearmesh stiffness. EJMC-A/Solids, 2009, 28(3): 461-468.

  c_f = cos²(αF) / (E·B) · [ L*(uf/Sf)² + M*(uf/Sf) + P*(1 + Q*·tan²αF) ]

几何量 (Chaari 2009 Fig.2 的标准定义):
  uf = 接触点到齿根圆的径向距离 (rc − rf)
  Sf = 齿根圆弦齿厚 = 2·rf·sin(θf)
  θf = 齿根半角 (取基圆半齿角 φb, 因基圆以下齿廓径向延伸)
  hfi = rf/r_int (齿根半径/轮缘内径)
L*,M*,P*,Q* 为 Sainsot et al. (2004) 多项式系数(θf 与 hfi 的函数)。
"""
import numpy as np

# Sainsot/Chaari 表1: 行 = [L*, M*, P*, Q*], 列 = [A B C D E F]
_COEF = np.array([
    [-5.574e-5,  -1.9986e-3, -2.3015e-4,  4.7702e-3,  0.0271, 6.8045],
    [ 60.111e-5,  28.100e-3, -83.431e-4, -9.9256e-3,  0.1624, 0.9086],
    [-50.952e-5, 185.50e-3,   0.0538e-4, 53.300e-3,  0.2895, 0.9236],
    [ -6.2042e-5,  9.0889e-3, -4.0964e-4,  7.8297e-3, -0.1472, 0.6904],
])


def chaari_foundation_compliance(g, r_contact, alphaF, r_int, B, E):
    theta_f = float(g.phi_b)                  # 齿根半角
    rf = g.rf
    if not (np.isfinite(r_int) and 0 < r_int < rf):
        r_int = rf / 2.0                      # 实体近似兜底(已在参数中登记[D]假设)
    hfi = rf / r_int
    uf = max(r_contact - rf, 0.0)
    Sf = 2 * rf * np.sin(theta_f)

    A, Bc, C, D, Ee, F = _COEF.T
    stars = A / theta_f**2 + Bc * hfi**2 + C * hfi / theta_f \
        + D / theta_f + Ee * hfi + F
    Ls, Ms, Ps, Qs = stars

    ratio = uf / Sf
    bracket = Ls * ratio**2 + Ms * ratio + Ps * (1 + Qs * np.tan(alphaF)**2)
    c_f = np.cos(alphaF)**2 / (B * E) * bracket

    detail = dict(theta_f=theta_f, rf=rf, r_int=r_int, hfi=hfi, uf=uf, Sf=Sf,
                  Lstar=Ls, Mstar=Ms, Pstar=Ps, Qstar=Qs, bracket=bracket,
                  c_f=c_f)
    return float(c_f), detail
