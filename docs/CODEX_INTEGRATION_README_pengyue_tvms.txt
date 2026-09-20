彭悦第3章 TVMS 集成说明（物理 SP + 物理 PR）
================================================

本路线只处理“势能法 TVMS -> 18DOF 动力学”链路。旧版
profileFloor + profileSin2*sin(...)^2 仅保留给 legacy 路线；当
p.mesh.tvmsSource='pengyue' 时，SP 和 PR 都使用预计算的物理单齿对
柔度并按活动齿对直接并联，不做二次正弦调制、插值平滑或频谱注入。

新增文件
--------
pg_build_pengyue_sp_tvms_model.m
  外啮合太阳轮-行星轮势能法 lookup；裂纹只修改指定故障齿的弯曲/剪切。
pg_eval_pengyue_sp_tvms.m
  按 epsilon 和 Parker gamma_sp^i=Zs*phi_p^i 评估 k_sp1..k_sp3。
pg_build_pengyue_pr_tvms_model.m
  内啮合行星轮-齿圈 lookup；行星侧沿用外齿势能积分，齿圈侧实现彭悦
  式(2-25)~(2-27)的 kb/ka/ks 内齿积分。
pg_eval_pengyue_pr_tvms.m
  按 gamma_rp^i=Zr*phi_p^i 评估 k_rp1..k_rp3，修复旧 ringOffset
  只改齿号而不改 pairAge 的相位错误。

构造与接入
----------
    p = pg_parameters_pengyue_fig314;
    fault = struct('type','none','tooth',1,'targetPlanet',1, ...
                   'crackDepth',0,'crackAngle',0);
    [p,spModel] = pg_prepare_pengyue_tvms(p,fault,'FoundationModel','none');
    [kSp,dSp] = pg_eval_pengyue_sp_tvms(t,p,fault,spModel);
    [kRp,dRp] = pg_eval_pengyue_pr_tvms(t,p,fault,p.mesh.pengyuePrModel);

pg_tvms 在 pengyue 模式下直接把 kSp/kRp 传给原 18DOF 装配；质量矩阵、
支承、坐标和动力学方程保持不变。

基础柔度边界
------------
彭悦表2-1没有给出太阳轮/行星轮的 r_int、theta_f，也没有给出齿圈
轮缘/轮体、齿根圆角几何。因此正式路线使用 FoundationModel='none'，
对应 c_f=0（kf=Inf），并在 model.ringFoundationStatus/report 中写明
MISSING。Chaari/Sainsot 模型只有在用户提供明确几何后才能启用；Chaari
的 u_f 定义也必须显式选择。旧的 7e-9 m/N 常数和 GlobalScale 校准值
已从正式路线禁用，传入非零常数会直接报错。

内啮合几何边界
--------------
论文未列出 alpha_f、theta、phi_b 的数值。当前代码从报告的 base/
addendum/root 半径导出可复现的标准渐开线近似，并将
model.internalGeometry.geometryStatus 标为 inferred；这不是 CAD 级
几何复现，正式论文应在获得齿圈齿根/轮缘图纸后替换这些字段。

健康验证
--------
运行 scripts/validate_pengyue_full_tvms.m 会输出 DOF、SP/PR 接触比、
三行星相位、单/双齿比例、相位回归误差和六条 TVMS CSV。对于
zs=17,zp=37,zr=91，理论相位为 SP=[0,2/3,1/3]，
PR=[0,1/3,2/3] cycles。

在健康 TVMS 通过几何审查之前，不应运行裂纹频谱，也不应调整 FFT、边带
或裂纹削弱幅度。
