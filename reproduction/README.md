# 彭悦第3章对照复现代码

本目录保存两套互相配套、但与仓库主模型隔离的复现代码：

- `pengyue_chapter3_reproduction/`：Python 势能法 TVMS、18DOF 动力学和
  数据导出；
- `my_gearbox_figures/`：读取导出 CSV 的 MATLAB 论文绘图脚本。

## 运行数据生成

从仓库根目录执行：

```powershell
python -m reproduction.pengyue_chapter3_reproduction.validation.run_user_gearbox_B
```

该入口使用用户齿轮参数 21/31/84、600 r/min，并将 `time_*.csv`、
`spec_*.csv` 和 TVMS 检查数据写入：

```text
reproduction/my_gearbox_figures/data/
```

`spec_*.csv` 是对应完整原始 DTE 响应去均值后的直接单边 FFT。生成脚本
不会在信号中人工增加啮合谐波、边带或冲击。

## 运行 MATLAB 出图

在 MATLAB 中执行：

```matlab
cd('reproduction/my_gearbox_figures')
addpath('scripts')
run_all_figures
```

生成的 PNG、TIF、PDF、EPS 和 FIG 位于该目录的 `export_*` 和
`fig_matlab/` 子目录。所有结果目录均已加入 `.gitignore`。

## 时域与频域定义

- `*_raw_paired` 时域图：各工况原始 DTE 响应，仅减去各自均值；
- 频谱图：同一原始 DTE 响应的单边幅值谱；
- `fault_component`：故障响应减健康响应的诊断残差，不作为原始系统响应。

健康/最小参数使用灰色，参数增加依次使用浅蓝、深蓝、橙色和红色。颜色只
表示工况顺序，不改变曲线数据。

## 证据边界

Python 对照实现用于检查TVMS、故障周期和频谱生成机制。它不替代仓库主线
MATLAB模型，也不包含齿圈模态传递路径或传感器观察器。论文引用结果前应
记录具体提交版本、运行参数和生成CSV的校验信息。
