# 无转速计多通道—多谐波联合相位解耦原型

## 方法主线

双通道第 `h` 个啮合谐波的解析信号写成

`z_s,h(t) = A_s,h(t) exp{i[h phi_m(t) + psi_s,h(theta_c(t))]}`，

其中 `phi_m` 是所有通道和谐波共享的啮合相位，`theta_c = phi_m/Zr`
是行星架相位，`psi_s,h` 是移动传递路径相位。对于三个等间隔行星轮，路径旋量采用
三行星空间阶次约束：

`P_s,h(theta) = sum_k c_s,h,k exp(i 3 k theta)`。

算法从双通道振动的多个固定频带出发，先以多通道、多谐波相位增量估计公共啮合相位，
再以通道互相位直接拟合差模路径旋量，并在零均值规范下把相对路径相位对称分配给两路。
路径补偿只乘单位模旋量，不修改实测幅值。机械零位的改变只会旋转傅里叶系数，因此算法
不需要机械零位、行星轮编号或故障齿编号。

## 可辨识性边界

绝对常相位不可由振动自身确定。双通道互相位能够识别的是差模移动路径；对所有通道完全
相同的公共路径相位仍位于公共相位规范中，不能在无额外观测时与转速相位无条件区分。
当前原型利用转速相位的跨通道、跨谐波共模性和路径的三行星空间周期性稳定估计。论文中
必须明确这一边界，并通过同步扭振或强载荷波动工况测试失效范围。

## 运行

```matlab
cd('<path-to>/planetary-gearbox-dynamic-model')
setup_paths
report = run_tacholess_joint_phase_demo
```

仿真真值和量化编码器计数只保存在 `report.truth` 中用于事后评分，不会传入估计器。主要
输出为：

- `results/tacholess_joint_phase_metrics.csv`
- `results/tacholess_joint_phase_validation.png`
- `results/tacholess_joint_phase_waveforms.png`
- `results/tacholess_joint_phase_demo.mat`

MAT 文件同时保存按估计行星架相位重采样的原始/校正周期矩阵，可直接接到 TSA、TPSVD
或知识角域包络流程。
