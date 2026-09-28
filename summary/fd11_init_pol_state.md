# 显式设置初始化后的偏振状态（2026-09-29）

`src/dev/fd11.cpp` 新增 `set_init_pol_state()`，供调用者在光子初始化后、第一次推进轨道前调用。
`main()` 中有两个调用位置：

```cpp
pe.init_random_radial(energy, initial_mode);
pe.set_init_pol_state();
```

```cpp
pe.perform_scattering(); // 内部调用 init()，建立新的光子路径
pe.set_init_pol_state();
```

`init()` 仍负责复位为 Mode、清空上一段偏振信息并初始化测地线 stepper。
stepper 初始化时已计算并缓存 FD11 事件的符号。
`set_init_pol_state()` 在偏振启用、当前为 Mode、缓存符号非负时立即启动积分；
否则保留初始化后的 Mode 状态。重复调用不会重新初始化已启动的偏振状态。

论文判据仍为 $\ell_A/r>\epsilon$，其中 $\ell_A=1/\kappa$。
乘以正数 $\kappa r$ 后等价于 $1-\epsilon\kappa r>0$；
代码仍将等号作为启动边界，阈值由 `--couple` 设置，默认为 `1e-3`。
本次只调整初始化流程，没有新增物理公式或更改判据。

为复用状态设置逻辑，`start_polarization(path, state, initial_h)` 接收明确的路径位置、
轨道状态和偏振初始步长。初始化调用使用 `x_old`、`y_old` 和 `h_old`，
不读取尚未准备的 dense output、`y_new` 或 `h_new`。
运行中事件仍通过无参 `start_polarization()` 使用已定位的步末位置与建议步长。
两种入口均设置本征模 Stokes、Integrating 阶段、启动位置、步长和计数，
并把启动事件的缓存符号设置为 -1，防止关闭事件时产生额外过零。

`step_geodesic()` 现在仅执行：

```cpp
stepper.do_step(0.1 * stepper.y_old[0]);
return stepper.detect_event();
```

删除了零长度启动事件以及相应的 `y_new`、`dydx_new`、`h_old` 人为重置。
外部若直接使用 `PhotonEvolution::init()` 或 `perform_scattering()`，也应在初始化之后
调用 `set_init_pol_state()`，再执行轨道演化。所有实际推进轨道的回归测试均同步使用这一流程。

验证命令：`make fd11-test`。测试覆盖初始未满足、初始已满足、重复调用、禁用偏振、
散射复位、后续事件定位以及原有独立 Jones 解/冻结/光学深度回归。
初始已满足的检查还验证：轨道位置和初始步长不变、本征模 Stokes 正确、事件符号已关闭。
日志保存在 `output/fd11_init_pol_tests.log`。

另运行 1000 光子、种子 4193 的完整输运检查，输出与日志分别为
`output/fd11_init_pol_smoke.txt`、`output/fd11_init_pol_smoke.log`。

验证结果：回归测试全部通过；完整输运逃逸 992、吸收 8，输出有限值和计数检查通过，最大 Stokes 范数平方误差 4.707e-13。`git diff --check` 通过。
