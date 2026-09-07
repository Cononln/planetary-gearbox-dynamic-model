# 配置与本地数据

`event_alignment_v1_plan.json` 是早期事件对齐路线的历史计划，包含本机数据路径、
旧路线参数和探索性阈值。它不是当前太阳轮带符号波形入口的配置文件，也不代表
里面每一项计划都已实现。保留它用于追溯，不直接当成新实验的已冻结协议。

需要运行该历史路线时，先复制成 `event_alignment.local.json`，再修改 `records`
中的路径并核对入口实际读取的参数。`*.local.json` 已被 Git 忽略。数据留在本机，
复制的配置不得放入密钥或账号令牌。

当前太阳轮入口用命令行传路径，示例：

```powershell
python scripts/test_sun_signed_waveform.py --source "D:/your-data/Sun050S600L00.tdms.mat"
```

运行前在仓库根目录创建 `results`，或先执行 `python scripts/run_smoke_tests.py`。
历史 MATLAB 数据入口也有本机默认路径；逐个查看函数参数，不假定所有入口接受
同一种配置。MATLAB 模型核心只依赖已提交的 `data/ring_modes_30.csv`，不需要原始测量。
