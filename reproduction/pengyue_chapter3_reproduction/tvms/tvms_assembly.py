# -*- coding: utf-8 -*-
"""
TVMS 装配: 单/双齿对并联 (式(2-24)/(3-3)) + Parker 相位 + 故障齿时序。

啮合运动学 (齿圈固定):
  太阳轮相对行星架转频 f_sm − f_c, 啮合频率 f_m = (f_sm−f_c)·zs
  行星轮 i 的局部啮合计数 u_i(t) = f_m·t + γ_sp,i/(2π)·zs ... 以"齿对序号"计:
    tooth-pair 序号 n = floor(u_i), 齿对内进度 ξ = u_i − n ∈ [0,1)
  每个齿对的实际啮合进度 η = ξ·ε (一个啮合周期内该齿对只在 ξ∈[0,ε) 承载,
  其余时间已脱开; 前一齿对在 ξ∈[0,ε−1) 仍在接触 → 双齿区)。

太阳轮裂纹齿 (编号 fault_tooth, 从1起): 当齿对序号 n ≡ fault_tooth−1 (mod zs)
时该齿对采用裂纹刚度。太阳轮相对行星架每转一周, 裂纹齿依次与三个行星轮
啮合一次 → 三次故障事件, 间隔 Ts = 1/(Np·f_m/zs) = 0.02967s。
"""
import numpy as np
from ..parameters import params as P


class TVMSModel:
    """预计算 k_pair(η) 查找表 (健康/裂纹), 供 ODE 与后处理快速插值。"""

    def __init__(self, comp_table_healthy, comp_table_fault, contact_ratio,
                 n_lookup=2001, fault_tooth=1):
        self.eta = np.linspace(0, 1, n_lookup)
        self.k_healthy = 1.0 / np.interp(self.eta, comp_table_healthy['eta'],
                                         comp_table_healthy['c_total'])
        if comp_table_fault is not None:
            self.k_fault = 1.0 / np.interp(self.eta, comp_table_fault['eta'],
                                           comp_table_fault['c_total'])
        else:
            self.k_fault = self.k_healthy.copy()
        self.epsilon = contact_ratio
        self.fault_tooth = fault_tooth
        self.n_lookup = n_lookup

    def k_pair(self, eta_deg, healthy_only=False):
        """齿对进度 η∈[0,1] 的单对齿刚度, 含故障齿判断。"""
        eta = np.clip(eta_deg, 0.0, 1.0)
        if healthy_only:
            return np.interp(eta, self.eta, self.k_healthy)
        return np.where(eta_deg >= 0,
                        np.interp(eta, self.eta, self.k_fault),
                        np.interp(eta, self.eta, self.k_healthy))

    def mesh_stiffness(self, u, healthy_only=False):
        """局部啮合计数 u (单位: 齿对周期)处的综合啮合刚度 (单/双齿对并联)。"""
        u = np.asarray(u, dtype=float)
        n = np.floor(u)
        xi = u - n
        k = np.zeros_like(u)
        active = []      # (n_pair, eta)
        # 当前齿对: ξ∈[0,ε) 承载, 进度 η=ξ/ε
        # 前一齿对: ξ∈[0,ε−1) 仍承载, 进度 η=(1+ξ)/ε
        e = self.epsilon
        m1 = xi < e
        eta1 = xi[m1] / e
        n1 = n[m1]
        # 齿对序号 n 对应太阳轮齿 = mod(n, zs)+1 (SP 啮合)
        if healthy_only:
            k[m1] = np.interp(eta1, self.eta, self.k_healthy)
        else:
            fault1 = np.mod(n1, P.zs) + 1 == self.fault_tooth
            kk = np.interp(eta1, self.eta, self.k_healthy)
            kk[fault1] = np.interp(eta1[fault1], self.eta, self.k_fault)
            k[m1] = kk
        active.append((n1, eta1, np.mod(n1, P.zs) + 1))
        if e > 1:
            m2 = xi < e - 1
            eta2 = (1.0 + xi[m2]) / e
            n2 = n[m2] - 1
            if healthy_only:
                k[m2] += np.interp(eta2, self.eta, self.k_healthy)
            else:
                fault2 = np.mod(n2, P.zs) + 1 == self.fault_tooth
                kk = np.interp(eta2, self.eta, self.k_healthy)
                kk[fault2] = np.interp(eta2[fault2], self.eta, self.k_fault)
                k[m2] += kk
            active.append((n2, eta2, np.mod(n2, P.zs) + 1))
        return k


def local_mesh_counter(t, planet_index, mesh='sp'):
    """行星轮 i 的局部啮合计数 (齿对周期数, 未取模)。

    γ_sp,i = +zs·φ_i (行星轮系内顺时针自转约定, 式(2-29));
    相位换算为齿对序偏移: Δu_i = γ_sp,i / (2π) (单位: 啮合周期)。
    """
    if mesh == 'sp':
        du = P.zs * P.phi_p[planet_index] / (2 * np.pi)
    else:
        du = -P.zr * P.phi_p[planet_index] / (2 * np.pi)
    return P.f_m * t + du


def sun_fault_event_times(planet_index, fault_tooth=None, n_events=6):
    """裂纹齿(太阳轮第 fault_tooth 齿)开始进入行星轮 i 啮合的时刻序列。

    本实现约定齿对啮合起点与局部计数整数对齐 (u_i = f_m·t + a_i, a_i = zs·φ_i/2π),
    故齿对 n 承载裂纹齿 ⟺ mod(n, zs) = fault_tooth−1。
    事件间隔 = zs/f_m (太阳轮相对行星架一周), 三行星轮事件互相错开 zs/(3·f_m)。
    """
    if fault_tooth is None:
        fault_tooth = P.fault_tooth
    a_i = P.zs * P.phi_p[planet_index] / (2 * np.pi)
    n_first = np.mod(fault_tooth - 1 - a_i, P.zs)      # 首个事件对应的 u−a 值 ∈[0,zs)
    t0 = n_first / P.f_m
    return t0 + P.zs / P.f_m * np.arange(n_events)
