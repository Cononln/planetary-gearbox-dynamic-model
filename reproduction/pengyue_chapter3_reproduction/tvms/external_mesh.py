# -*- coding: utf-8 -*-
"""
势能法单对齿啮合柔度 —— 外啮合 (太阳轮-行星轮), 论文式(2-6)~(2-23)。

直接按应变能积分实现 (不使用论文式(2-18)~(2-20)的α参数化闭式, 二者数学等价,
积分形式可避免文本提取造成的符号歧义):
  1/kb = ∫ [cosαF·(ξF−ξ) − sinαF·hF]² / (E·I(ξ)) dξ      式(2-6)
  1/ka = ∫ sin²αF / (E·A(ξ)) dξ                            式(2-7)
  1/ks = ∫ 1.2·cos²αF / (G·A(ξ)) dξ                        式(2-8)
  kh   = πEB / (4(1−ν²))                                    式(2-21)
  kf   = Chaari/Sainsot 齿基柔度, 式(2-22)                  文献[71]
单对齿总柔度 (式(2-23)):
  1/k = 1/kh + Σ_gear (1/kb + 1/ka + 1/ks + 1/kf)

裂纹 (式(3-1)(3-2)): 只修改指定齿轮的弯曲/剪切截面:
  hx ≤ hq:  Ix = (2hx)³B/12,  Ax = 2hx·B            (健全)
  hx > hq:  Ixc = (hx + hc − q·sinγ)³B/12,
            Axc = (hx + hc − q·sinγ)·B              (有效全齿厚削减)
  hq = hc − q·sinγ, hc = 齿根半齿厚
赫兹/轴向压缩/齿基/配对齿轮保持健全 (论文3.2.1节明确)。
"""
import numpy as np
from ..parameters import params as P
from .geometry import GearGeometry
from .foundation import chaari_foundation_compliance


def sun_geometry():
    return GearGeometry(P.zs, P.rs, P.rb_s, P.ra_s, P.rf_s, P.alpha0)


def planet_geometry():
    return GearGeometry(P.zp, P.rp, P.rb_p, P.ra_p, P.rf_p, P.alpha0)


def ring_geometry():
    return GearGeometry(P.zr, P.rr, P.rb_r, P.ra_r, P.rf_r, P.alpha0, internal=True)


def tooth_compliance(g: GearGeometry, r_contact, r_int, crack=None, n_int=400):
    """单个齿轮轮齿在接触半径 r_contact 处的弯曲/轴向/剪切柔度 (每单位啮合力)。

    crack = None 或 dict(q=裂纹深度m, gamma=裂纹角rad)
    返回 (c_bending, c_axial, c_shear, detail)
    """
    rb, rf, alpha0 = g.rb, g.rf, g.alpha0
    alphaF = np.arccos(min(max(rb / r_contact, 0.0), 1.0))
    hF = float(g.half_thickness(r_contact))     # 接触点到中心线距离 (式2-11)
    xiF = float(g.xi(r_contact))                # 接触点中心线坐标 (式2-12)

    # 沿中心线积分: 从齿根截面到接触截面
    r = np.linspace(rf, r_contact, n_int)
    hx = g.half_thickness(r)
    xi = g.xi(r)

    s_full = 2.0 * hx                          # 健全全齿厚
    if crack is not None and crack.get('q', 0) > 0:
        q = crack['q']
        gamma = crack['gamma']
        hc = g.half_thickness_at_root()
        hq = hc - q * np.sin(gamma)
        s_crack = hx + hc - q * np.sin(gamma)  # 式(3-1)(3-2)有效全齿厚
        s_eff = np.where(hx > hq, s_crack, s_full)
        s_eff = np.minimum(s_eff, s_full)      # 数值保护(条件已保证)
        s_eff = np.maximum(s_eff, 1e-6)        # 防退化
    else:
        hc = hq = np.nan
        s_eff = s_full

    B, E, G = P.B_face, P.E, P.G
    I = B * s_eff**3 / 12.0
    A = B * s_eff

    # 弯矩 M/F = cosαF·(ξF−ξ) − sinαF·hF   (式2-6, 力臂为截面到接触点中心线距离)
    M_over_F = np.cos(alphaF) * (xiF - xi) - np.sin(alphaF) * hF

    # 去重排序保证 trapz 单调
    order = np.argsort(xi)
    xi_s, M_s = xi[order], M_over_F[order]
    I_s, A_s = I[order], A[order]

    c_bend = float(np.trapezoid(M_s**2 / (E * I_s), xi_s))
    c_axial = float(np.trapezoid(np.sin(alphaF)**2 / (E * A_s), xi_s))
    c_shear = float(np.trapezoid(1.2 * np.cos(alphaF)**2 / (G * A_s), xi_s))

    # Chaari 齿基柔度 (式2-22, 文献[71])
    c_f, f_detail = chaari_foundation_compliance(
        g, r_contact, alphaF, r_int, P.B_face, P.E)

    detail = dict(alphaF=alphaF, hF=hF, xiF=xiF, hc=hc, hq=hq,
                  s_root=float(s_eff[0]), s_root_healthy=float(s_full[0]),
                  foundation=f_detail)
    return c_bend, c_axial, c_shear, c_f, detail


def external_pair_compliance(eta, crack_sun=None, crack_planet=None,
                             r_int_sun=None, r_int_planet=None, n_int=400,
                             alpha_op=None, a_op=None):
    """太阳轮-行星轮单对齿柔度随啮合进程 η∈[0,1] (η=0 啮入, 1 啮出)。

    alpha_op/a_op: 工作压力角与工作中心距(宽安装/变位时≠标准值),
    默认取 params 的标准值。齿厚几何仍按参考渐开线(x=0)。
    """
    if r_int_sun is None:
        r_int_sun = P.r_int_sun
    if r_int_planet is None:
        r_int_planet = P.r_int_planet
    alpha_op = P.alpha0 if alpha_op is None else alpha_op
    a_op = (P.rs + P.rp) if a_op is None else a_op

    gs, gp = sun_geometry(), planet_geometry()
    a = a_op
    # 啮合线几何: l_s + l_p = a'·sinα' (工作压力角)
    l_sum = a * np.sin(alpha_op)
    l_s_tip = np.sqrt(gs.ra**2 - gs.rb**2)     # 太阳轮齿顶极限
    l_p_tip = np.sqrt(gp.ra**2 - gp.rb**2)     # 行星轮齿顶极限
    # 啮入(η=0): 行星轮齿顶与太阳轮齿根侧接触 → l_s = l_sum − l_p_tip
    # 啮出(η=1): 太阳轮齿顶与行星轮齿根侧接触 → l_s = l_s_tip
    l_s_start = l_sum - l_p_tip
    l_s_end = l_s_tip

    eta = np.atleast_1d(np.asarray(eta, dtype=float))
    l_s = l_s_start + eta * (l_s_end - l_s_start)
    l_p = l_sum - l_s

    r_cs = np.sqrt(gs.rb**2 + l_s**2)          # 太阳轮上接触点半径
    r_cp = np.sqrt(gp.rb**2 + l_p**2)

    rows = []
    for i in range(eta.size):
        cs = tooth_compliance(gs, r_cs[i], r_int_sun, crack=crack_sun, n_int=n_int)
        cp = tooth_compliance(gp, r_cp[i], r_int_planet, crack=crack_planet, n_int=n_int)
        c_h = 4 * (1 - P.nu**2) / (np.pi * P.E * P.B_face)   # 1/kh 式(2-21)
        c_tot = c_h + (cs[0] + cs[1] + cs[2] + cs[3]) + (cp[0] + cp[1] + cp[2] + cp[3])
        rows.append((eta[i], c_h, *cs[:3], cs[3], *cp[:3], cp[3], c_tot))
    cols = ['eta', 'c_hertz', 'sun_c_bend', 'sun_c_axial', 'sun_c_shear', 'sun_c_found',
            'pla_c_bend', 'pla_c_axial', 'pla_c_shear', 'pla_c_found', 'c_total']
    return {c: np.array([r[i] for r in rows]) for i, c in enumerate(cols)}
