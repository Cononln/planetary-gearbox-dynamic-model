# Planetary gearbox dynamic model

## 当前入口 / Start here

本仓库用于行星齿轮箱动力学、时变传递路径和多测点故障脉冲对齐研究。
这是研究代码，数值自检通过与实验方法有效是两件事。

当前主动力学入口是 **18 个集中参数自由度 + 7 个有效齿圈模态坐标**，
`paper_v3` 主入口的测点为 0°、120°、240°，敏感方向暂按径向向外。
独立的 ANSYS 齿圈后处理路径可以使用 0°、90°、120°、240° 四测点；这四个
FE 标签不等同于 18DOF 主入口的三通道配置。太阳轮输入、齿圈固定、行星架输出；
当前配置输入扭矩为 20 N m，对应理想功率平衡的行星架负载 100 N m。
模型路径、未标定参数和当前算法状态见 [CURRENT_STATE](docs/CURRENT_STATE.md)。

## 代码组织 / Code organization

仓库保持一个轻量、可追溯的研究代码布局：

- `src/`：可复用的 MATLAB 模型函数（参数、矩阵、TVMS、动力学和传递路径）；
- `scripts/`：按研究步骤运行的入口、数据准备和诊断脚本；
- `src_py/`：可复用的 Python 信号/相位算法与 FE 资产读取函数；
- `tests/`：不依赖个人数据的 MATLAB/Python 合同与数值自检；
- `configs/`、`docs/`、`reproduction/`：参数说明、方法边界和隔离复现路线；
- `results/`、`data/`：本地运行目录，生成结果和原始数据不作为源码模块提交。

MATLAB 入口先调用根目录的 `setup_paths`，Python 入口从 `scripts/` 调用
`src_py/` 的函数。新增的齿圈模态传递路径也遵循这套分层：计算内核在
`src/`，运行、诊断图和文件导出在 `scripts/`；每个入口的输入/输出和适用假设写在
对应的 `docs/` 说明中。

FE 路径的坐标重注册、逐事件齿序和模态贡献诊断均为只读后处理工具：它们读取
已保存的 18DOF 源响应和用户提供的 ANSYS 资产，不修改源动力学求解。原始 FE
文件和生成结果保留在本地，跨机器运行时应通过命令行显式传入路径。

### MATLAB 模型

在 MATLAB 中把当前目录切换到仓库根目录，然后执行：

```matlab
setup_paths
report = validate_paper_v3_model
```

需要生成响应时再执行（不是上传或安装时自动运行的步骤）：

```matlab
output = simulate_paper_v3_three_channel_response
```

默认输出健康、太阳轮 25% 和 50% 断齿三种工况，采样率 51200 Hz，
保留 3 s，另计算并丢弃前 0.2 s。结果写入本地 `results/`。
该入口使用支持 `arguments` 块的现代 MATLAB；全部历史脚本可能还需要
Signal Processing Toolbox。本次仓库整理不重新运行完整动力学实验。

### Python 算法与快速自检

建议使用 Python 3.11 或更新版本，在仓库根目录执行：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe scripts/run_smoke_tests.py
```

快速自检只使用已知合成信号，检查 Pearson 时移符号、SVD/TSA 等价关系、
带符号波形的谱幅值与抵消、太阳轮/行星轮潮汐周期定义，不读取实测数据。
它不验证物理模型标定或实验诊断优越性。测试文件使用独立脚本入口，
不要把普通 `pytest` 的收集结果当成这些自检已经执行。

太阳轮实测入口需要本地数据，显式传入路径：

```powershell
.\.venv\Scripts\python.exe scripts/test_sun_signed_waveform.py --source "D:/your-data/Sun050S600L00.tdms.mat"
```

该脚本按既有记录格式读取 `chanvals`：前三列振动、后两列编码器，采样率
按既有记录约定使用 51200 Hz，分析第 30–90 s；这不是通用数据导入器。
它会加载完整 MAT 数组，需要足够内存。其他记录格式先阅读入口代码和
[配置说明](configs/README.md)，不能仅修改文件名就假定列映射一致。

### 多 AI 协作与文件范围

- 每项修改使用独立分支和独立克隆目录，合并前看差异并运行相关测试。
- 新 AI 先读 [AGENTS.md](AGENTS.md)、[当前状态](docs/CURRENT_STATE.md)
  和 [协作流程](docs/COLLABORATION.md)。
- 仓库保留源码、模型参数、技术说明和必要的 30 阶模态频率 CSV。
  原始采集数据、CAD、论文全文、生成图表和本地运行环境保留在本机。
- 私有仓库并不会自动授权其他 AI 服务访问；连接服务或邀请协作者由用户决定。
- 当前未授予开源许可证，研究材料的对外发布须由所有者决定。

## 独立的 ANSYS 齿圈模态路径

18DOF 集中参数模型产生三路带符号 PR 动态啮合力；新增的有限元路径模块
将它们投影到 ANSYS 导出的齿圈模态，再恢复固定测点径向响应。
该模块不把齿圈模态坐标并入 18DOF 源方程，也不覆盖上述 18+7 有效模态
研究入口。代码分工、模态资产准备、运行命令和证据边界见
[齿圈有限元模态传递路径](docs/fe_ring_modal_transfer_path.md)。
原始 ANSYS 文件和生成结果不随 Git 仓库上传。当前 CSV→SI 转换后的模态
质量归一化尚待从 ANSYS 求解设置或广义质量核实，因此该路径的绝对幅值
暂作探索性结果，不作为已标定的 m/s² 定量结论。

## 历史工作流 / Legacy reference

下文保留早期双通道和方法原型的说明，便于追溯；其中的 “default” 指对应
历史版本。当前三通道研究以顶部 v3 入口和 CURRENT_STATE 为准。

This MATLAB project builds the planetary-gearbox model used to test the
dual-channel phase-alignment method.  The default response generator is now
a two-dimensional 18-DOF lumped-parameter model (LPM).  The former elastic
ring/modal workflow is retained only as an optional comparison.

## Default lumped-parameter workflow

The sun, ring, carrier, and three planets each have `[x,y,theta]`, giving 18
gearbox DOFs.  Tooth-level TVMS and transmission error generate the dynamic
mesh forces.  The physical contact forces recovered from the LPM drive two
local SDOF mass-spring-damper path observers at 0 and 90 degrees through
smooth moving-source proximity windows.  No FE mode shapes or elastic-ring
modal coordinates enter the default response simulation.

Run the default model from the project directory:

    setup_paths
    reports = run_all_validations
    output = simulate_sensor_responses
    lowFrequency = plot_low_frequency_spectrum

The complete equations and observation model are documented in
`docs/default_lpm_sensor_model.md`.

## Loaded v2 comparison workflow

The provisional publication-oriented comparison model adds balanced
drive/load torques and replaces the local SDOF observers with the optional
six-mode elastic-ring coordinates.  It is kept separate from the original
baseline so that assumptions remain auditable:

    setup_paths
    output = simulate_loaded_v2_responses
    report = validate_loaded_v2_response

See `docs/loaded_v2_model.md`.  The present 100 N m carrier load is an
explicit placeholder, not a measured operating condition.

## Stages

1. 18-DOF lumped-parameter gear model.
2. Partial-tooth-break time-varying mesh stiffness (TVMS).
3. Frequency-calibrated analytical ring modes.
4. Time-varying paths to the fixed 0 and 90 degree sensors.
5. Phase-difference and group-delay truth fields.

Run the optional legacy elastic-ring validation from the project directory:

    setup_paths
    reports = run_hybrid_validations

All SI-unit parameters and assumptions are defined in src/pg_parameters.m.
The bearing and mean mesh stiffness values are initial engineering estimates
and are explicitly marked as calibration parameters.

## Optional root-crack fault model

A Shen-style single-sun-tooth root crack is available alongside the
partial-width broken-tooth baseline. It computes an energy-based reduction
of the damaged tooth's bending and shear stiffness for crack extensions
0.5, 1.5, and 2.5 mm while leaving Hertz and axial terms unchanged:

    setup_paths
    report = validate_root_crack_model
    output = simulate_root_crack_responses(20)

See `docs/root_crack_fault_model.md`. The reference paper does not report
the imposed speed-ripple amplitude, so that optional term defaults to zero.

## Coordinate convention

The bodies are ordered as sun, ring, carrier, planet 1, planet 2, planet 3.
Each body has [x, y, theta], so q has 18 entries. Translation is in m and
rotation in rad. The dynamic coordinates are perturbations about the known
steady motion. The mean carrier angle is prescribed as
phi_c(t) = phi_c0 + omega_c*t, while carrier theta remains a torsional
perturbation DOF.

The ring is not removed from the model: its three DOFs are retained and its
mean rotation is fixed through the housing support. This is required before
the elastic ring response is coupled in Stage 3.

## Evidence boundary

The supplied modal-frequency text contains frequencies but no FE mode-shape
vectors. Stage 3 will therefore use a frequency-calibrated analytical ring
basis, not claim an FE modal reduction.

## Main outputs

- results/all_lpm_validation_reports.mat: default 18-DOF LPM checks.
- results/lpm_phase_landmarks.csv: 0/90-degree amplitudes and cross-phases at
  the mesh, path-sideband, and fault-related frequencies.
- results/carrier_order_sidebands.png/csv: `fm+n*fc` carrier-modulation
  components and their 0/90-degree cross-phases.
- results/mesh_band_carrier_envelope.png/csv: demodulated 2/4/6-Hz carrier
  orders from the 150-186-Hz band around the mesh fundamental.
- results/step2_tvms.png: healthy/25%/50% tooth-level TVMS.
- results/step3_ring_modes.png: analytical mode shapes and calibrated FRF.
- results/step4_paths.png: moving-path magnitude and dual-sensor phase.
- results/phase_truth_field.mat: optional hybrid-model complex FRF truth.
- results/phase_truth_band_delay.csv: carrier-angle delay lookup table.
- results/step5_phase_truth.png: truth-field diagnostic figure.
- results/sensor_response_time.png: healthy/25%/50% dual-sensor waveforms.
- results/sensor_response_spectrum.png: corresponding spectra to 5.5 kHz.
- results/sensor_response_sun_cases.mat/csv: numerical acceleration signals.
- results/sensor_response_time_1s.png: one-second dual-sensor waveforms.
- results/sensor_response_spectrum_1s.png: one-second-record spectra.
- results/sensor_response_spectrum_log_1s.png: log-scale diagnostic spectra.
- results/sensor_response_envelope_spectrum_1s.png: 24 Hz fault modulation.
- results/sensor_response_spectrum_0_200Hz.png: mesh frequency and sidebands.
- results/sensor_response_sun_cases_1s.mat/csv: one-second numerical signals.

The optional hybrid-model group delay is an effective structural group delay
that includes resonance and anti-resonance energy storage. It is not a pure
material time-of-flight.

## Tacholess joint phase prototype

The blind dual-channel/multiharmonic prototype separates common speed phase
from moving-path phase with three-planet spatial orders and applies
phase-only compensation:

    setup_paths
    report = run_tacholess_joint_phase_demo

Encoder phase is hidden from the estimator and retained only for scoring.
See `docs/tacholess_joint_phase_method.md` for the model, identifiability
boundary, baselines, and TSA/TPSVD interface.

The measured 50%-sun-tooth-break TDMS file can be processed with both
tachometer columns hidden until scoring:

    report = run_real_tacholess_joint_phase

The runner reads vibration columns 3/4, estimates the mesh center from
vibration alone, and audits the three-planet spatial-order assumption using
the encoder only after the blind result has been frozen.

Measured-data findings and the resulting soft-3k/asymmetry refinement are
documented in `docs/real_tacholess_phase_results.md`.

To inspect the known 50% sun-tooth-break feature after phase alignment:

    diagnosis = run_real_phase_aligned_diagnosis

This compares single-channel, naive dual-channel, and phase-aligned fusion
in the envelope spectrum and carrier-cycle TPSVD domains.

## Three-channel measured-data prototype

The legacy MAT records contain three synchronized vibration channels at
0, 120, and 240 degrees plus carrier-encoder A/B phases.  The complete
pairwise phase-graph prototype separates common-phase sensor selection from
differential-path edge reliability:

    report = run_real_three_channel_phase( ...
        fullfile('data_local','sun','Sun050S600L00.tdms.mat'))

Low-coherence path edges are excluded without altering their measured
amplitudes.  Encoder samples remain hidden until scoring.  Initial P050 and
Sun050 findings, including the present failure to recover an encoder-angle
path oracle, are documented in `docs/real_three_channel_phase_results.md`.

The measured-data runner now adds an alternating-carrier-cycle safeguard.
It validates no correction, every two-sensor edge, and the complete graph on
held-out carrier revolutions.  A complete graph that reduces validation
coherence automatically falls back to the most coherent stable pair or to
no correction.  The selector ranks absolute held-out coherence after a
positive-gain admission test, so a very weak edge is not rewarded merely for
starting near zero.

For the 2026 BL600/PF50 records, run the common/differential diagnostic:

    diagnosis = run_three_channel_safe_diagnosis

The aligned common mode must not be used as the only diagnostic signal.
The phase-aligned common response and the differential residual are retained
together for phase-synchronous TPSVD.  In the current BL/PF50 comparison the
common-only fault contrast decreases, while the differential residual fault
family increases by 2.78 dB; this is the present evidence for the
common/differential interface rather than direct equal-weight averaging.

## 彭悦（2025）图 3-14 专项复现入口

为避免污染当前 21/31/84 三通道研究模型，新增了独立的彭悦论文图 3-14 复现配置：

```matlab
setup_paths
report = validate_pengyue_fig3_14
output = reproduce_pengyue_fig3_14
```

该入口使用论文的 17/37/91 齿数、800 r/min、50 N·m、20 μm 传递误差、
q=4 mm 和 γ=15/30/45/60/75°，并输出太阳轮—行星轮 DTE 的五层三维频谱。
论文未报告的齿宽、位移无量纲标尺、绝对 TVMS 等不会被伪装成已知值；当前显式
假设和第一版尚未加入的 backlash/逐齿摩擦项见
[专项复现说明](docs/pengyue_fig3_14_reproduction.md)。

当前用于用户 21/31/84 齿轮参数的 Python 18DOF 对照复现和 MATLAB 论文
绘图源码位于 [reproduction](reproduction/README.md)。该目录只跟踪模型、
验证和绘图代码；生成的 CSV、PNG、FIG、PDF 等结果由 `.gitignore` 排除。
