# -*- coding: utf-8 -*-
"""
内啮合 (行星轮-内齿圈) 势能法 TVMS —— 论文式(2-25)~(2-27)。

内啮合几何: 两基圆内切, 啮合线上切点间距
  T_rT_p = a·sinα0,  且 l_r = l_p + a·sinα0 (接触点在两切点外侧延长段)
接触区: l_p 从 (l_r,ringTip − a·sinα0) 到 行星轮齿顶极限。
行星轮为主动轮(齿圈固定、行星架驱动内啮合): 啮入时接触点靠近行星轮齿根,
与论文3.2.3节"内啮合过程行星轮作为主动轮, 啮入阶段靠近故障齿齿根"一致。

内齿圈轮齿: 齿厚随半径增大而增大 (θ(r) = π/(2z) − invα0 + invα(r)),
悬臂梁从齿根(外侧 rf)向齿顶(内侧 ra)方向承载。
"""
import numpy as np
from ..parameters import params as P
from .external_mesh import tooth_compliance, planet_geometry, ring_geometry


def internal_pair_compliance(eta, crack_planet=None, r_int_planet=None,
                             r_int_ring=None, n_int=400):
    """行星轮-内齿圈单对齿柔度随 η∈[0,1] (η=0 啮入)。

    内啮合啮入 = 行星轮齿根侧/内齿圈齿顶侧进入接触
    (行星轮主动, 与论文图3-6规律: 啮入阶段裂纹影响小)。
    """
    if r_int_planet is None:
        r_int_planet = P.r_int_planet
    if r_int_ring is None:
        r_int_ring = P.r_int_ring

    gp, gr = planet_geometry(), ring_geometry()
    a = P.rr - P.rp
    d_tan = a * np.sin(P.alpha0)              # 两切点沿啮合线间距
    l_p_tip = np.sqrt(gp.ra**2 - gp.rb**2)    # 行星轮齿顶极限
    l_r_tip = np.sqrt(gr.ra**2 - gr.rb**2)    # 内齿圈齿顶极限(内侧)
    # 啮入(η=0): 内齿圈齿顶与行星轮齿根侧接触 → l_p = l_r_tip − d_tan
    # 啮出(η=1): 行星轮齿顶与内齿圈齿根侧接触 → l_p = l_p_tip
    l_p_start = l_r_tip - d_tan
    l_p_end = l_p_tip

    eta = np.atleast_1d(np.asarray(eta, dtype=float))
    l_p = l_p_start + eta * (l_p_end - l_p_start)
    l_r = l_p + d_tan

    r_cp = np.sqrt(gp.rb**2 + l_p**2)
    r_cr = np.sqrt(gr.rb**2 + l_r**2)

    rows = []
    for i in range(eta.size):
        cp = tooth_compliance(gp, r_cp[i], r_int_planet, crack=crack_planet, n_int=n_int)
        cr = tooth_compliance(gr, r_cr[i], r_int_ring, crack=None, n_int=n_int)
        c_h = 4 * (1 - P.nu**2) / (np.pi * P.E * P.B_face)
        c_tot = c_h + (cp[0] + cp[1] + cp[2] + cp[3]) + (cr[0] + cr[1] + cr[2] + cr[3])
        rows.append((eta[i], c_h, *cp[:3], cp[3], *cr[:3], cr[3], c_tot))
    cols = ['eta', 'c_hertz', 'pla_c_bend', 'pla_c_axial', 'pla_c_shear', 'pla_c_found',
            'ring_c_bend', 'ring_c_axial', 'ring_c_shear', 'ring_c_found', 'c_total']
    return {c: np.array([r[i] for r in rows]) for i, c in enumerate(cols)}
