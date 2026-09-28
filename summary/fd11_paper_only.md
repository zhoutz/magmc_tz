# 仅保留 FD11 论文启动判据（2026-09-29）

按用户要求，从 `src/dev/fd11.cpp` 删除方向转动率启动条件、专用 `chi_prime()`、
`adiabatic_event()`、判据枚举、选项字段和 `--coupling-criterion` 命令行分支。
现在直接运行 `build/fd11` 即使用论文判据；`--couple` 仍可调整阈值，默认 `1e-3`。

论文式 (25) 给出 $\ell_A=1/\kappa$，故对于 $\kappa>0,r>0$，

\[
\ell_A/r>\epsilon
\iff 1/(\kappa r)>\epsilon
\iff 1>\epsilon\kappa r
\iff 1-\epsilon\kappa r>0.
\]

`coupling_event()` 只计算最后一个表达式；当 $\kappa=0$ 时也能返回有限值 1。
仍由 geodesic stepper 定位零点，在 `evolve_geodesic()` 末尾启动偏振积分。
初始已满足条件、散射后重新检测和关闭偏振的处理继续保留。

Stokes 积分器中随磁场转动坐标基所需的角度导数属于演化方程，
不是此次删除的启动判据，因而继续保留。冻结和光学深度算法未改动。

同步清理了测试中的双判据接口及方向导数启动测试；其余轨迹对照均使用 FD11 判据。
原双判据性能报告保留为历史记录，当前使用说明和公式推导首页已更新。

验证：

- `make fd11-test` 通过，包括启动事件、散射复位、独立 Jones 解、冻结及光学深度测试。
  日志：`output/fd11_paper_only_tests.log`。
- 使用删除前源码编译 `build/fd11_before_remove_adiabatic`，选择 `--coupling-criterion fd11`，
  与当前程序对照：1000 光子、种子 4193、默认能量与容差。
  两者均逃逸 992、吸收 8、散射 971 次；能量输出完全相同，方向余弦最大差为 2.165e-15。
  输出不逐字节相同，Stokes 向量最大差为 1.2025e-5；接受步数为 73013 与 73006，
  拒绝步数均为 8199。此次代码简化未改变 FD11 的数学表达式，但不承诺自适应浮点轨迹逐位相同。
  独立 Jones 回归测试同时通过；这里的差值是两份默认容差计算的差，而非对真实解误差的严格上界。
- 当前样本最大 $|Q^2+U^2+V^2-1|=2.518e-13$，输出有限值检查通过。
- 已移除的 `--coupling-criterion` 返回 Unknown option。
- `git diff --check` 通过。

原始对照输出及日志：`output/fd11_paper_only_before.*`、`output/fd11_paper_only_after.*`；
差值记录：`output/fd11_paper_only_comparison.json`。本次未重新进行性能基准测试。
