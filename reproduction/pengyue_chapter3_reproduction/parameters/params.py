# -*- coding: utf-8 -*-
"""
彭悦(2025)《行星齿轮故障动力学模型与特征研究》第3章复现 —— 参数定义。

参数分类标注:
  [A] 论文表2-1/正文明确给出
  [B] 可由论文参数计算
  [C] 论文引用文献可恢复
  [D] 论文未给出, 采用标准工程假设(必须在 unknown_parameters.md 中登记)
  [E] 当前完全无法确定

所有单位均为 SI (m, kg, N, s, rad), 除非另行注明。
"""
import numpy as np

# ----------------------------------------------------------------齿轮几何 [A]/[B]
zs, zp, zr = 17, 37, 91            # [A] 表2-1 齿数
Np = 3                              # [A] 行星轮个数
m_mod = 5e-3                        # [A] 表2-1 模数 m/mm = 5
alpha0 = np.deg2rad(20.0)           # [A] 标准压力角20° (用户指令确认; 论文未单列, 属[A]级共识)

rs = m_mod * zs / 2                 # [B] 太阳轮分度圆半径 42.5 mm
rp = m_mod * zp / 2                 # [B] 行星轮分度圆半径 92.5 mm
rr = m_mod * zr / 2                 # [B] 内齿圈分度圆半径 227.5 mm

# 表2-1 基圆直径 [A] (mm): 79.8793 / 173.8431 / 427.5601
rb_s = 79.8793e-3 / 2               # = rs*cos(20°) 验证一致
rb_p = 173.8431e-3 / 2
rb_r = 427.5601e-3 / 2

# 齿顶/齿根圆 [D]: 论文只给 m 与 z, 未给齿顶高系数; 采用标准齿制
#   外齿轮: ha*=1.0, c*=0.25 → ra = r + m, rf = r - 1.25m
#   内齿圈: 齿向内, ra = r - m, rf = r + 1.25m
ha_star, c_star = 1.0, 0.25
ra_s = rs + ha_star * m_mod
rf_s = rs - (ha_star + c_star) * m_mod
ra_p = rp + ha_star * m_mod
rf_p = rp - (ha_star + c_star) * m_mod
ra_r = rr - ha_star * m_mod         # 内齿圈齿顶(朝内)
rf_r = rr + (ha_star + c_star) * m_mod

rc_carrier = rs + rp                # [B] 行星轮中心半径 135 mm
d_bc_carrier = 253.7e-3             # [A] 行星架(等效)基圆直径, 用于扭转支撑换算

# ----------------------------------------------------------------材料 [D]
E  = 206e9                          # [D] 弹性模量(钢), 论文未给出
nu = 0.30                           # [D] 泊松比, 论文未给出
G  = E / (2 * (1 + nu))             # [B] 剪切模量

# 齿宽 [D]: 论文未给出。势能法所有柔度项均 ∝1/B, 故 k_pair ∝ B 严格成立,
# B 同时承担绝对标定。默认值由 validation/fig3_2 的敏感性表锁定后写回。
B_face = 6.0e-3                     # 锁定值(见 results/face_width_sensitivity.csv)

# Chaari/Sainsot 齿基几何 [D]: 论文未给轮缘内径 r_int。
# 假设轮缘厚度 rim = rf - r_int, 对齿基多项式影响弱(见敏感性表)。
rim_sun    = 10 * m_mod             # 太阳轮轮缘厚 50 mm → r_int_s = -13.75 mm<0 → 用实体近似
rim_planet = 10 * m_mod
rim_ring   = 10 * m_mod
# r_int 必须满足 0 < r_int < rf; 太阳轮 rf=36.25mm, 取 r_int = rf/2 (轴孔量级) [D]
r_int_sun    = max(rf_s - rim_sun,    rf_s / 2)
r_int_planet = max(rf_p - rim_planet, rf_p / 2)
r_int_ring   = max(rf_r - rim_ring,   rf_r / 2)

# ----------------------------------------------------------------质量/惯量 [A] 表2-1
ms, mp_, mr, mc = 20.29, 3.91, 59.88, 57.314          # kg
Is = 27049.57e-6      # kg·m²  (论文单位 kg·mm² → ×1e-6)
Ip = 22368.33e-6
Ir = 3853845.79e-6
Ic = 776910.239e-6

# ----------------------------------------------------------------支撑 [A] 表2-1
kx_s, kx_p, kx_r, kx_c = 10e6, 58e6, 1000e6, 195e6    # N/μm → N/m
kt_s_lin, kt_r_lin, kt_c_lin = 10e6, 57.6e6, 10e6     # N/μm(基圆切向线刚度约定)
zeta_support = 0.03                                   # [A] 支撑阻尼系数
zeta_mesh    = 0.07                                   # [A] 啮合阻尼比

# ----------------------------------------------------------------工况 [A] 3.4.1
ns_rpm   = 800.0                    # [A] 太阳轮输入转速 rpm
Tout     = 50.0                     # [A] 输出负载 N·m (用户指令; 论文3.4.1文字为5500N·mm, 登记)
e_av     = 20e-6                    # [A] 综合啮合误差幅值 m
b_half_dimless = 3.0                # [A] 无量纲半齿侧间隙 (物理换算需 b_c, 论文未给 → [E])
mu_fric  = 0.06                     # [A] 摩擦系数

# 裂纹基准工况 [A] 图3-2/3-9/3-10
q_crack   = 3e-3                    # m
gamma_deg = 45.0

# 故障齿编号 [D]: 论文未指明裂纹齿相位, 取 fault_tooth=15 使三行星轮裂纹事件
# 在图3-3的16周期窗口内的位置(≈14.0/8.33/2.67 周期)与论文布局一致。
# 该选择只平移事件时刻, 不影响 fs 谱结构与事件间隔。
fault_tooth = 15

# ----------------------------------------------------------------特征频率 [B] 表3-1
f_sm = ns_rpm / 60.0                                   # 13.333 Hz 太阳轮转频
# 行星架: ωc = ωs·zs/(zs+zr) (齿圈固定)
f_c  = f_sm * zs / (zs + zr)                           # 2.099 Hz  (表3-1: 2.1 Hz)
f_m  = (f_sm - f_c) * zs                               # 191.11 Hz (表3-1: 191.1 Hz)
f_s_fault = f_m * Np / zs                              # 33.72 Hz  (表3-1: 33.7 Hz)
f_p_fault = f_m / zp                                   # 5.166 Hz  (表3-1: 5.2 Hz)
T_s = 1.0 / f_s_fault                                  # 0.02967 s (论文: 0.0297 s)

# 行星轮安装角 [B]: 等间距 0°,120°,240° (论文未显式给, 由 Np=3 及图2-1)
phi_p = 2 * np.pi * np.arange(Np) / Np

# Parker 啮合相位 [C] 式(2-29)/(2-30), 行星轮在行星架系内顺时针自转约定:
#   γ_sp_i = +zs·φ_i  (mod zs),  γ_rp_i = −zr·φ_i  (mod zr)
gamma_sp_cycles = np.mod(zs * phi_p / (2 * np.pi), 1.0)   # 啮合周期分数相位
gamma_rp_cycles = np.mod(-zr * phi_p / (2 * np.pi), 1.0)


def summary():
    s = []
    s.append(f"epsilon_SP = {contact_ratio('sp'):.6f}")
    s.append(f"epsilon_PR = {contact_ratio('pr'):.6f}")
    s.append(f"f_m = {f_m:.4f} Hz, f_s = {f_s_fault:.4f} Hz, f_p = {f_p_fault:.4f} Hz, f_c = {f_c:.4f} Hz")
    s.append(f"T_s = {T_s:.5f} s (论文 0.0297 s)")
    s.append(f"Parker SP phase (cycles) = {gamma_sp_cycles}")
    s.append(f"Parker RP phase (cycles) = {gamma_rp_cycles}")
    return "\n".join(s)


def contact_ratio(which):
    pb = np.pi * m_mod * np.cos(alpha0)   # 基节
    if which == 'sp':
        a = rs + rp
        path = (np.sqrt(ra_p**2 - rb_p**2) + np.sqrt(ra_s**2 - rb_s**2)
                - a * np.sin(alpha0))
    else:
        a = rr - rp
        path = (np.sqrt(ra_p**2 - rb_p**2) - np.sqrt(ra_r**2 - rb_r**2)
                + a * np.sin(alpha0))
    return path / pb


if __name__ == "__main__":
    print(summary())
