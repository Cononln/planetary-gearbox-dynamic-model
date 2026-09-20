# -*- coding: utf-8 -*-
"""
彭悦第2章 18 自由度平移-扭转耦合行星齿轮动力学模型 (式2-37~2-46)。

自由度 (行星架旋转坐标系 OXY; 低速忽略陀螺/离心项, 式(2-42)):
  q = [xs,ys,θs, xr,yr,θr, xc,yc,θc, x1,y1,θ1, x2,y2,θ2, x3,y3,θ3]

啮合投影 (式2-1/2-3):
  δ_sp_i = u_i·(q_s − q_pi) + rbs·θs − rbp·θpi − e_sp_i,  u_i = (cos(φi+α), sin(φi+α))
  δ_pr_i = v_i·(q_pi − q_r) + rbp·θpi − rbr·θr − e_pr_i,  v_i = (cos(φi−α), sin(φi−α))
啮合力 (式2-33): F = k(t)·f(δ) + c·δ̇, f 为齿侧间隙非线性 (式2-32)
摩擦 (式2-34~2-36): F_f = μ·sgn(R)·F 作用于垂直啮合线方向,
  力臂 = 接触点到齿轮中心的啮合线坐标 (符号在节点翻转)。
"""
import numpy as np
from ..parameters import params as P

DOF = 18


class Model18:
    def __init__(self, k_sp_lookup, k_pr_lookup, epsilon_sp, epsilon_pr,
                 a_sp, a_rp, fault_sp_tooth=None, ksp_cracked=None,
                 fault_pr=None, backlash_m=0.0, friction=True,
                 e_av=None, sp_fault_gear='sun', alpha_sp=None, alpha_pr=None):
        self.eps_sp, self.eps_pr = epsilon_sp, epsilon_pr
        self.eta_sp, self.ksp_h = k_sp_lookup[0], k_sp_lookup[1]
        self.ksp_c = ksp_cracked
        self.eta_pr, self.kpr_h = k_pr_lookup[0], k_pr_lookup[1]
        self.a_sp = np.asarray(a_sp, float)
        self.a_rp = np.asarray(a_rp, float)
        self.fault_tooth = fault_sp_tooth
        self.sp_fault_gear = sp_fault_gear  # 'sun': 齿号按zs计; 'planet': 按zp计
        self.fault_pr = fault_pr            # dict(planet, tooth, k, offset_cycles=18.5)
        self.backlash = backlash_m
        self.friction = friction
        self.e_av = P.e_av if e_av is None else e_av

        m = np.zeros(DOF)
        m[0] = m[1] = P.ms; m[2] = P.Is
        m[3] = m[4] = P.mr; m[5] = P.Ir
        m[6] = m[7] = P.mc; m[8] = P.Ic
        for i in range(3):
            m[9 + 3 * i] = m[10 + 3 * i] = P.mp_; m[11 + 3 * i] = P.Ip
        self.M = m

        K = np.zeros(DOF)
        K[0] = K[1] = P.kx_s; K[2] = P.kt_s_lin * P.rb_s ** 2
        K[3] = K[4] = P.kx_r; K[5] = P.kt_r_lin * P.rb_r ** 2
        K[6] = K[7] = P.kx_c; K[8] = P.kt_c_lin * (P.d_bc_carrier / 2) ** 2
        for i in range(3):
            K[9 + 3 * i] = K[10 + 3 * i] = P.kx_p
        C = np.zeros(DOF)
        C[0] = C[1] = 2 * P.zeta_support * np.sqrt(K[0] * P.ms)
        C[2] = 2 * P.zeta_support * np.sqrt(K[2] * P.Is)
        C[3] = C[4] = 2 * P.zeta_support * np.sqrt(K[3] * P.mr)
        C[5] = 2 * P.zeta_support * np.sqrt(K[5] * P.Ir)
        C[6] = C[7] = 2 * P.zeta_support * np.sqrt(K[6] * P.mc)
        C[8] = 2 * P.zeta_support * np.sqrt(K[8] * P.Ic)
        for i in range(3):
            C[9 + 3 * i] = C[10 + 3 * i] = 2 * P.zeta_support * np.sqrt(P.kx_p * P.mp_)
        self.Ks, self.Cs = K, C

        m_eq_sp = 1.0 / (1.0 / P.ms + 1.0 / P.mp_)
        m_eq_pr = 1.0 / (1.0 / P.mp_ + 1.0 / P.mr)
        self.k_mean_sp = float(np.mean(self.ksp_h))
        self.k_mean_pr = float(np.mean(self.kpr_h))
        self.c_sp = 2 * P.zeta_mesh * np.sqrt(self.k_mean_sp * m_eq_sp)
        self.c_pr = 2 * P.zeta_mesh * np.sqrt(self.k_mean_pr * m_eq_pr)

        self.phi = P.phi_p.copy()
        a_sp_ = P.alpha0 if alpha_sp is None else alpha_sp
        a_pr_ = P.alpha0 if alpha_pr is None else alpha_pr
        self.u_sp = np.stack([np.cos(self.phi + a_sp_), np.sin(self.phi + a_sp_)])
        self.u_pr = np.stack([np.cos(self.phi - a_pr_), np.sin(self.phi - a_pr_)])
        self.rbs, self.rbp, self.rbr, self.rc = P.rb_s, P.rb_p, P.rb_r, P.rc_carrier

        self.Tin = P.Tout * P.zs / (P.zs + P.zr)   # 理想功率平衡 (D类假设, 登记)
        self.Tout = P.Tout
        self.mu = P.mu_fric

        # 啮合线几何 (从基圆切点量起)
        self.lsun_pitch = P.rs * np.sin(P.alpha0)
        self.lsun_start = np.sqrt(P.ra_s ** 2 - P.rb_s ** 2) - (
            np.sqrt(P.ra_p ** 2 - P.rb_p ** 2) - self.lsun_pitch)
        self.lsun_span = np.sqrt(P.ra_s ** 2 - P.rb_s ** 2) - self.lsun_start
        self.lpla_pitch = P.rp * np.sin(P.alpha0)
        self.lpla_start = np.sqrt(P.ra_r ** 2 - P.rb_r ** 2) - (
            P.rr - P.rp) * np.sin(P.alpha0)
        self.lpla_span = np.sqrt(P.ra_p ** 2 - P.rb_p ** 2) - self.lpla_start

        self.te_ph_sp = 2 * np.pi * np.mod(P.zs * self.phi / (2 * np.pi), 1.0)
        self.te_ph_pr = 2 * np.pi * np.mod(-P.zr * self.phi / (2 * np.pi), 1.0)

    # ------------------------------------------------刚度
    def k_sp(self, t, i):
        u = P.f_m * t + self.a_sp[i]
        n = np.floor(u); xi = u - n
        e = self.eps_sp
        z_fault = P.zs if self.sp_fault_gear == 'sun' else P.zp
        k = np.interp(xi / e, self.eta_sp, self.ksp_h)
        if self.ksp_c is not None:
            k = np.where(np.mod(n, z_fault) + 1 == self.fault_tooth,
                         np.interp(xi / e, self.eta_sp, self.ksp_c), k)
        if xi < e - 1:
            eta_b = (1 + xi) / e
            kb = np.interp(eta_b, self.eta_sp, self.ksp_h)
            if self.ksp_c is not None:
                kb = np.where(np.mod(n - 1, z_fault) + 1 == self.fault_tooth,
                              np.interp(eta_b, self.eta_sp, self.ksp_c), kb)
            k = k + kb
        return k

    def k_pr(self, t, i):
        u = P.f_m * t + self.a_rp[i]
        n = np.floor(u); xi = u - n
        e = self.eps_pr
        use_fault = (self.fault_pr is not None and i == self.fault_pr['planet'])
        off = float(self.fault_pr.get('offset_cycles', P.zp / 2.0)) if use_fault else 0.0
        nf = np.floor(u + off)          # fault-tooth cycle index, offset-aware
        k = np.interp(xi / e, self.eta_pr, self.kpr_h)
        if use_fault:
            k = np.where(np.mod(nf, P.zp) + 1 == self.fault_pr['tooth'],
                         np.interp(xi / e, self.eta_pr, self.fault_pr['k']), k)
        if xi < e - 1:
            eta_b = (1 + xi) / e
            kb = np.interp(eta_b, self.eta_pr, self.kpr_h)
            if use_fault:
                kb = np.where(np.mod(nf - 1, P.zp) + 1 == self.fault_pr['tooth'],
                              np.interp(eta_b, self.eta_pr, self.fault_pr['k']), kb)
            k = k + kb
        return k

    # ------------------------------------------------摩擦臂 (接触点啮合线坐标)
    def _arm(self, u, e, l_start, l_span):
        xi = u - np.floor(u)
        eta = np.where(xi < e, xi / e, 0.0)
        return l_start + eta * l_span

    # ------------------------------------------------右端项
    def rhs(self, t, q, v):
        F = -self.Ks * q - self.Cs * v
        F[2] += self.Tin
        F[8] -= self.Tout

        xc, yc, thc = q[6], q[7], q[8]
        vc, vyc, vthc = v[6], v[7], v[8]
        for i in range(3):
            cphi, sphi = np.cos(self.phi[i]), np.sin(self.phi[i])
            xci = xc - thc * self.rc * sphi
            yci = yc + thc * self.rc * cphi
            vxci = vc + vthc * self.rc * sphi
            vyci = vyc - vthc * self.rc * cphi
            kp, cp = P.kx_p, 2 * P.zeta_support * np.sqrt(P.kx_p * P.mp_)
            fx = -kp * (q[9 + 3 * i] - xci) - cp * (v[9 + 3 * i] - vxci)
            fy = -kp * (q[10 + 3 * i] - yci) - cp * (v[10 + 3 * i] - vyci)
            F[9 + 3 * i] += fx; F[10 + 3 * i] += fy
            F[6] -= fx; F[7] -= fy
            F[8] += fx * self.rc * sphi - fy * self.rc * cphi

        for i in range(3):
            # ---- SP ----
            us = self.u_sp[:, i]
            w = 2 * np.pi * P.f_m
            d = (us[0] * (q[0] - q[9 + 3 * i]) + us[1] * (q[1] - q[10 + 3 * i])
                 + self.rbs * q[2] - self.rbp * q[11 + 3 * i]
                 - self.e_av * np.sin(w * t + self.te_ph_sp[i]))
            dd = (us[0] * (v[0] - v[9 + 3 * i]) + us[1] * (v[1] - v[10 + 3 * i])
                  + self.rbs * v[2] - self.rbp * v[11 + 3 * i]
                  - self.e_av * w * np.cos(w * t + self.te_ph_sp[i]))
            fn = _gap(d, self.backlash)
            Fm = float(self.k_sp(t, i)) * fn + self.c_sp * dd
            F[0] -= Fm * us[0]; F[1] -= Fm * us[1]; F[2] -= Fm * self.rbs
            F[9 + 3 * i] += Fm * us[0]; F[10 + 3 * i] += Fm * us[1]
            F[11 + 3 * i] += Fm * self.rbp
            if self.friction:
                l_s = float(self._arm(P.f_m * t + self.a_sp[i], self.eps_sp,
                                      self.lsun_start, self.lsun_span))
                l_p = float(self._arm(P.f_m * t + self.a_sp[i], self.eps_sp,
                                      self.lsun_start, self.lsun_span))  # 同一接触点
                R = l_s - self.lsun_pitch
                Ff = self.mu * abs(Fm) * np.sign(R)
                tx, ty = -us[1], us[0]
                F[0] += Ff * tx; F[1] += Ff * ty
                F[2] += -Ff * l_s
                F[9 + 3 * i] -= Ff * tx; F[10 + 3 * i] -= Ff * ty
                F[11 + 3 * i] += -Ff * l_p

            # ---- PR ----
            vp = self.u_pr[:, i]
            d = (vp[0] * (q[9 + 3 * i] - q[3]) + vp[1] * (q[10 + 3 * i] - q[4])
                 + self.rbp * q[11 + 3 * i] - self.rbr * q[5]
                 - self.e_av * np.sin(w * t + self.te_ph_pr[i]))
            dd = (vp[0] * (v[9 + 3 * i] - v[3]) + vp[1] * (v[10 + 3 * i] - v[4])
                  + self.rbp * v[11 + 3 * i] - self.rbr * v[5]
                  - self.e_av * w * np.cos(w * t + self.te_ph_pr[i]))
            fn = _gap(d, self.backlash)
            Fm = float(self.k_pr(t, i)) * fn + self.c_pr * dd
            F[9 + 3 * i] -= Fm * vp[0]; F[10 + 3 * i] -= Fm * vp[1]
            F[11 + 3 * i] -= Fm * self.rbp
            F[3] += Fm * vp[0]; F[4] += Fm * vp[1]; F[5] += Fm * self.rbr
            if self.friction:
                l_p2 = float(self._arm(P.f_m * t + self.a_rp[i], self.eps_pr,
                                       self.lpla_start, self.lpla_span))
                l_r = l_p2 + (P.rr - P.rp) * np.sin(P.alpha0)
                R = l_p2 - self.lpla_pitch
                Ff = self.mu * abs(Fm) * np.sign(R)
                tx, ty = -vp[1], vp[0]
                F[9 + 3 * i] += Ff * tx; F[10 + 3 * i] += Ff * ty
                F[11 + 3 * i] += Ff * l_p2
                F[3] -= Ff * tx; F[4] -= Ff * ty
                F[5] += -Ff * l_r
        return F / self.M

    # ------------------------------------------------DTE
    def dte_sp(self, t, q, i):
        us = self.u_sp[:, i]
        return (us[0] * (q[0] - q[9 + 3 * i]) + us[1] * (q[1] - q[10 + 3 * i])
                + self.rbs * q[2] - self.rbp * q[11 + 3 * i]
                - self.e_av * np.sin(2 * np.pi * P.f_m * t + self.te_ph_sp[i]))

    def dte_pr(self, t, q, i):
        vp = self.u_pr[:, i]
        return (vp[0] * (q[9 + 3 * i] - q[3]) + vp[1] * (q[10 + 3 * i] - q[4])
                + self.rbp * q[11 + 3 * i] - self.rbr * q[5]
                - self.e_av * np.sin(2 * np.pi * P.f_m * t + self.te_ph_pr[i]))


def _gap(d, b):
    if b <= 0:
        return d
    if d > b:
        return d - b
    if d < -b:
        return d + b
    return 0.0


def rk4_run(model, t_end, dt, q0=None, v0=None, record=True):
    n = int(round(t_end / dt))
    q = np.zeros(DOF) if q0 is None else q0.copy()
    v = np.zeros(DOF) if v0 is None else v0.copy()
    if record:
        T = np.zeros(n); Q = np.zeros((n, DOF))
    for j in range(n):
        t = j * dt
        if record:
            T[j] = t; Q[j] = q
        a1 = model.rhs(t, q, v)
        a2 = model.rhs(t + dt / 2, q + v * dt / 2, v + a1 * dt / 2)
        a3 = model.rhs(t + dt / 2, q + (v + a1 * dt / 2) * dt / 2,
                       v + a2 * dt / 2)
        a4 = model.rhs(t + dt, q + (v + a2 * dt / 2) * dt, v + a3 * dt)
        q = q + (v + (a1 + a2 + a3) * dt / 6) * dt
        v = v + (a1 + 2 * a2 + 2 * a3 + a4) * dt / 6
    if record:
        return T, Q, v
    return q, v
