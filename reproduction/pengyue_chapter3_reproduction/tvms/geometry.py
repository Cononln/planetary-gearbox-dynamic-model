# -*- coding: utf-8 -*-
"""
渐开线轮齿几何 (论文式(2-11)~(2-17) 的精确几何实现)。

坐标系: 齿中心线为 ξ 轴 (原点在齿轮中心), 接触点在齿槽法向截面内。
约定 (标准渐开线几何, 与论文式(2-17)一致):
  θ(r) = φb − inv(α(r))                r ≥ rb   (inv 为渐开线函数)
  hx(r) = r·sinθ(r)                    r ≥ rb   (半齿厚)
  基圆以下 (rf ≤ r < rb): 论文式(2-17)第一分支, hx = rb·sinφb 常数
  ξ(r) = r·cosθ(r)                              (沿中心线坐标)
其中 φb = π/(2z) + inv(α0) 为基圆半齿角。
"""
import numpy as np


def inv(a):
    """渐开线函数 inv(α) = tanα − α"""
    return np.tan(a) - a


class GearGeometry:
    def __init__(self, z, r_pitch, r_base, r_tip, r_root, alpha0, internal=False):
        self.z = z
        self.rp = r_pitch
        self.rb = r_base
        self.ra = r_tip
        self.rf = r_root
        self.alpha0 = alpha0
        self.internal = internal
        # 基圆半齿角 (外齿轮与内齿轮齿廓同为渐开线, 齿距角相同)
        self.phi_b = np.pi / (2 * z) + inv(alpha0)

    def theta(self, r):
        """半径 r 处的半齿角。基圆以下(外齿轮)取 φb(径向延伸, 对应式2-17分支1)。"""
        r = np.asarray(r, dtype=float)
        th = np.full_like(r, self.phi_b)
        mask = r > self.rb
        alpha_r = np.arccos(np.clip(self.rb / r[mask], -1, 1))
        th[mask] = self.phi_b - inv(alpha_r)
        # 内齿轮齿厚随半径增大而增大: θ(r) = π/(2z) − inv(α0) + inv(α(r))
        if self.internal:
            th = np.full_like(r, np.pi / (2 * self.z) - inv(self.alpha0))
            alpha_r = np.arccos(np.clip(self.rb / r, -1, 1))
            th = th + inv(alpha_r)
        return th

    def half_thickness(self, r):
        """半齿厚 hx(r)。外齿轮基圆以下按论文式(2-17)分支1取常数 rb·sinφb。"""
        r = np.asarray(r, dtype=float)
        if self.internal:
            return r * np.sin(self.theta(r))
        hx = np.where(r >= self.rb,
                      r * np.sin(self.theta(r)),
                      self.rb * np.sin(self.phi_b))
        return hx

    def xi(self, r):
        """沿齿中心线坐标 ξ(r) = r·cosθ(r)。"""
        return np.asarray(r, dtype=float) * np.cos(self.theta(r))

    def half_thickness_at_root(self):
        """裂纹起点 hc: 齿根表面到中心线距离 (式(3-1)的 hc)。"""
        return float(self.half_thickness(self.rf))
