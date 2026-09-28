# 齿圈有限元模态传递路径

本模块把 18 自由度集中参数动力学模型的行星轮—齿圈（PR）动态啮合力
作为移动法向载荷投影到 ANSYS 导出的齿圈模态，再恢复固定测点的径向响应。
动力学源模型与齿圈路径是两个明确的阶段：

    18DOF Newmark 响应
        → 带符号的三路 PR 啮合力和移动接触角
        → ANSYS 齿圈模态投影
        → 测点径向响应

## 代码分工

- `scripts/ansys_export_unified_fe_asset_dat2.txt`：从现有 ANSYS 结果集导出三测点 CSV；它不重新求解模态。
- `src_py/fe_ring_asset_io.py`：共用 CSV 读取与单位转换函数。
- `scripts/prepare_ansys_exact_three_sensor_asset.py`：0/120/240 三测点 CSV 到 MAT 的命令行入口。
- `scripts/prepare_ansys_exact_four_sensor_asset.py`：0/90/120/240 四测点 CSV 到 MAT 的命令行入口。
- src/pg_fe_ring_modal_model.m：读取质量归一化模态资产并检查字段和尺寸。
- src/pg_fe_ring_contact_angle.m：将源动力学角度注册到 FE 全局坐标。
- src/pg_fe_ring_input_shape.m：移动啮合位置的模态参与因子周期插值。
- src/pg_fe_ring_modal_path_gain.m：接触角到测点的复模态频响及 |H|² 诊断量。
- src/pg_apply_fe_ring_transfer_causal.m：逐时步精确零阶保持（ZOH）模态状态积分。
- src/pg_apply_fe_ring_transfer.m：保留的频域对照；正式因果事件图使用逐时步实现。
- scripts/run_causal_fe_ring_path_comparison.m：健康/太阳轮裂纹的因果路径对比入口。
- scripts/diagnose_fe_modal_path_gain.m：接触角和测点的模态参与诊断。

## 模态资产

运行顺序为：在 ANSYS 已求解的同一 RST 中执行 APDL 导出 → 选用对应
测点数量的 Python 构建器 → 执行 Python 轻量自检 → 在 MATLAB 中运行
因果路径入口。原始 RST、DAT、节点坐标和全量模态 CSV 不进入 Git。
三测点和四测点构建器是两种输入合同，按实际导出的测点数选择其一，
不是要求顺次运行两个命令：

    python scripts/prepare_ansys_exact_three_sensor_asset.py --source D:/path/to/ansys_export --output results/fe_ring_asset
    python scripts/prepare_ansys_exact_four_sensor_asset.py --source D:/path/to/ansys_export --output results/fe_ring_asset

仓库附带的 ANSYS APDL 导出示例在完整 FE 网格中寻找 0/120/240 测点，
并检查节点一一对应和不超过 2 mm 的最近节点映射误差。四测点构建器要求
本地导出文件还包含人工选定的 SENSOR_90 测点，且每阶四个测点的节点数一致；
该测点不能由三测点文件自动补造。
APDL 的 `auto_sensor_check.csv` 和完成标记用于检查导出状态；手工导出的
历史 CSV 可能没有这两个文件，此时构建器只做 CSV 内部一致性检查，
使用者仍需核对原始 RST 与节点选区。每次 APDL 重导出应使用独立目录，
避免把旧文件误当成同一次结果。

两个构建入口均接受 `mode_shapes_sensor.csv` 的两种列格式：

- `mode_id,sensor_id_or_angle,node_id,x_m,y_m,z_m,ux_m,uy_m,uz_m`（9 列，第二列可为 `SENSOR_120` 或数值 `120`）；
- `mode_id,frequency_Hz,angle_deg,node_id,x,y,z,ux,uy,uz`（10 列，既有 ANSYS 导出）。

当前 N-mm-MPa 导出中的坐标和模态位移会转换成 SI 米制；生成的 MAT 同时保留
`sensor_angle_local_deg` 和由节点坐标计算的 `sensor_angle_global_deg`。
`--ring-nodes` 和 `--sensor-nodes` 必须与实际导出的每阶节点数相符。
本地角度标签从 `SENSOR_0` 起按顺时针增加；FE 全局 `atan2` 角按逆时针增加。
构建器逐个检查标签和实际测点坐标，120°/240° 若交换会报错。旧版按逆时针
生成的 CSV 应保留原样，用相应真实角度重新核对标签后再构建，不可直接重贴标签。

四测点资产的本地标签可以是 0、90、120、240 度；资产同时保存 FE 全局角度，
路径计算通过 `pg_fe_ring_contact_angle` 注册源模型角度，不能把本地标签
直接当作全局坐标。MAT 文件被
.gitignore 排除，避免把个人有限元结果或生成数据误提交到仓库。
当前注册采用“源模型局部 0° 对应 `SENSOR_0` 方向”的坐标约定；
实际行星架装配零位和转向仍需与 FE 模型/试验记录核对，不能把约定当作
已测得的绝对相位。

## 运行示例

    setup_paths
    out = run_causal_fe_ring_path_comparison( ...
        'AssetPath', "results/fe_ring_asset/ansys_ring_modal_transfer_asset_exact_four_sensor.mat", ...
        'OutputDir', "results/fe_ring_path_demo", ...
        'RecordDuration', 2.0, ...
        'DiscardDuration', 0.2, ...
        'CrackDepthMm', 0.90);

该入口先求解 18DOF 响应，再把三个 PR 动态啮合力按瞬时接触位置投影到
全部保留模态。没有角距离经验增益、人工延迟、平滑或通道缩放。模态阻尼
和 ANSYS 模态归一化会影响绝对幅值。当前构建器把 N-mm-MPa 导出中的
振型分量按 mm→m 换算，并写入 `mass_normalized=1`，但尚未从 RST 广义质量
核实换算后是否仍为 SI 单位质量归一化；因此当前数值只能用于探索性比较，
不能直接作为已标定的 m/s² 或每牛顿传递率。应先核对 ANSYS 求解归一化及
质量单位、必要时修正模态尺度，再与实验或锤击 FRF 比较。
peak_lag_s 是事件窗内最大故障增量峰的位置，不是材料波飞行时间或首次到达时间。
若不指定 `AssetPath`，运行入口先查找四测点 MAT，找不到再查找三测点
MAT；概览图优先对比 0°/90°，三测点时对比 0°/120°。其余测点均保留在
CSV 和事件指标中。

## 检查

    setup_paths
    report = validate_paper_v3_model;
    smoke = test_fe_ring_modal_transfer;

    diagnostic = diagnose_fe_modal_path_gain( ...
        'AssetPath', "results/fe_ring_asset/ansys_ring_modal_transfer_asset_exact_four_sensor.mat");

从仓库根目录单独运行 Python 轻量自检：

    python scripts/run_smoke_tests.py

前一项检查验证 18DOF 运动学和矩阵维度，后一项使用合成小模型检查模态
方程、周期接触插值和零载荷响应；它们都不替代 ANSYS 资产的模态收敛性
或实验标定。Python 自检使用临时合成 CSV 核对列格式、测点角度和 MAT 字段，
不会启动完整动力学仿真。`diagnose_fe_modal_path_gain` 需已有真实 MAT 资产，
与上述无需资产的轻量自检不同。
