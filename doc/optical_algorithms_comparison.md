2026-09-27，本地工作区基于 `860b828`，将三个指定历史算法移植后，与 human、ref、test_qagp 在同一轮实验中重新比较。

**结论：这批样本中，默认综合选择是 global_sqrt；径向吞吐优先可选 transport_fast；非径向精度优先可选收紧容差的 velocity_panels。没有一种配置在所有指标上胜出。**

global_sqrt 默认相对于 human 的径向/随机非径向加速分别为 **2.47/1.86 倍**，相对于 ref 为 **8.08/8.58 倍**。它保持了接近 human 的径向精度，随机案例全部返回，但其中 4 个案例采用了历史实现的舍入警告处理策略，不能等同于全部通过严格 GSL 状态检查。transport_fast 的默认径向速度最快，精度稍低；在随机路径上出现求积失败，且比 global_sqrt 慢。

**移植内容与一致性**

| 新文件 | 原始提交 | 保留的核心机制 |
|---|---|---|
| `src/bench/global_sqrt.cpp` | `a907a14c98e03c8183a514791c2303a302d15bbe` | 缓存整条轨道、合并物理事件区间、平方根换元、稳定的 D/u² 计算、极值辅助搜索 |
| `src/bench/velocity_panels.cpp` | `9204013bc5919c2edad4e53c3c58bfe8cde0a4db` | 分布分位点与速度事件、D=0 分段、sin² 端点换元、8/16 点 Gauss 自适应求积 |
| `src/bench/transport_fast.cpp` | `4be9c249402dab313c171bc04244bc2be6a39563` | 历史 DOPR5 控制器、热分布速度标记、直接 GK(7,15) 求积、径向几何缓存及尾部裁剪 |

三份 `.cpp` 含算法实现和独立可执行入口。共同输入/输出位于 `benchmark_io.hpp`；transport_fast 独有的历史步进器和求积器保留在 `src/bench/detail/`，避免替换当前项目的同名组件。

适配仅涉及当前 `Boltzmann(b0,n_knots)` / `b_mean` 接口、当前位置的 Photon 初值及测试入口。global_sqrt/velocity_panels 使用当前 DOPR5 自动准备 dense output，因此去掉了旧代码中的重复准备调用。transport_fast 保留自己的历史步进器、默认容差、`tail_u=12`、`split_field_knots=false`、`guard_extremum_bounds=false`，没有把审计配置冒充历史默认配置。

将三个历史提交的原始数值代码分别提取到隔离目录，用相同 1900 个输入独立编译运行：**三个移植版与对应历史版的成功/失败状态完全相同，全部成功结果的最大绝对差均为 0**。transport_fast 的 3 个默认随机失败也在历史代码中复现。三个新程序的无参数径向输出与批量测试接口逐项一致；之前三个原程序的径向输出也再次验证一致。

**测试协议与参考解**

- 900 个径向案例：`beta0=-0.1,...,-0.9`，`muz=0,...,0.9`，能量为 `0.01,0.1,1,10,100`，每个组合测 E/O 两种偏振，从 `r0=10` 向外发射。
- 1000 个随机非径向案例：500 条独立路径，每条 E/O 两种偏振；种子 `20260927`。beta0 均匀分布于 `[-0.9,-0.1]`，muz 均匀分布于 `[-0.999,0.999]`，切向方位角均匀分布于 `[0,2π]`，能量对数均匀分布于 `[0.01,100]`。
- 前 250 条路径从表面发射，cos(alpha) 均匀分布于 `[0,1]`；后 250 条的 r0 对数均匀分布于 `[10.01,100]`，cos(alpha) 均匀分布于 `[-1,1]`。传播至吸收 `r=10` 或逃逸 `r=10000`。
- 两个输入文件与上一轮 `build/qagp_compare` 的文件逐字节相同。human/ref/test_qagp 使用上一轮相同的通用路径测试适配；特别是非径向 ref 不再于首次 D=0/beta=0 提前终止。
- 本机 macOS arm64、GCC 16.2.0、`-std=c++23 -O3`、Homebrew GSL。每种配置完整预热一次，之后串行交替顺序运行 5 次，取整批求解耗时中位数。六种方法均重新测量，没有拼接上一轮计时。
- 主表包含每条路径的分布构造/速度标记准备、轨道、事件和求积，排除编译、启动、磁场表加载与文件输出；含轻量计数。失败尝试的耗时保留，失败结果不计入误差。额外报告 velocity_panels 缓存分位点的模式。

径向参考仍为 `table/bench_od.txt`，并重新用独立速度空间积分核验：收紧容差后最大相对变化 `5.46e-13`，与该表最大相对差 `7.36e-13`。它不依赖被测算法的轨道积分或空间共振求根。

非径向参考为 a907a14 中的细步长 `ref_regular`，进行平方根换元。coarse/fine 的 `(最大步长/r, 轨道容差, 求积容差)` 分别为 `(0.00125,1.25e-12,1.25e-9)` 和 `(0.0003125,3.125e-13,3.125e-10)`。全部随机案例最大绝对变化 `1.63e-8`、最大相对变化 `2.08e-8`。再针对六种算法及对照配置的最坏结果、失败案例和参考变化选出 **63 个案例**，以 `(0.000078125,7.8125e-14,1e-12)` 加密，与 fine 最大绝对差 `1.51e-9`、最大相对差 `4.15e-10`。各参考档位最终无失败或求积警告。

非径向参考仍是数值收敛参照，且与 global_sqrt 共享若干物理、事件定位和换元代码，不是独立解析真值。这限制了对共同系统误差的排除能力。以下相对误差仅对 `|tau_ref|>1e-10` 计算；绝对误差包含全部成功结果。

**默认配置：精度、速度及返回状态**

| 样本 | 算法 | 总耗时（秒） | 成功/总数 | 最大绝对误差 | 最大相对误差 | P99 相对误差 |
|---|---|---:|---:|---:|---:|---:|
| 径向 | human | 0.1997 | 900/900 | 6.60e-8 | 5.48e-9 | 3.29e-9 |
| 径向 | ref | 0.6532 | 900/900 | 8.29e-9 | 4.71e-9 | 2.35e-9 |
| 径向 | test_qagp | 0.2661 | 900/900 | 6.60e-8 | 5.48e-9 | 3.29e-9 |
| 径向 | global_sqrt | 0.08079 | 900/900 | 6.54e-8 | 5.43e-9 | 3.28e-9 |
| 径向 | velocity_panels | 1.0330 | 900/900 | 6.53e-8 | 5.42e-9 | 3.28e-9 |
| 径向 | transport_fast | **0.05453** | 900/900 | 4.20e-6 | 3.48e-7 | 7.55e-8 |
| 随机非径向 | human | 0.4059 | 999/1000 | 2.33e-7 | 1.55e-6 | 3.82e-8 |
| 随机非径向 | ref | 1.8685 | 1000/1000 | 1.02e-6 | 1.59e-6 | 4.17e-8 |
| 随机非径向 | test_qagp | 0.5834 | 998/1000 | 2.31e-7 | 1.55e-6 | 3.79e-8 |
| 随机非径向 | global_sqrt | **0.2177** | 1000/1000* | 3.47e-7 | 2.24e-6 | 4.58e-8 |
| 随机非径向 | velocity_panels | 1.8315 | 1000/1000 | 4.44e-6 | 1.52e-6 | 1.01e-7 |
| 随机非径向 | transport_fast | 0.5340 | 997/1000 | 3.19e-6 | 9.03e-7 | 1.07e-7 |

*“成功”表示按各自历史错误处理规则返回有限非负结果，不代表使用了统一严格容差验收。global_sqrt 的警告接受规则详见下文。

默认全体径向结果表明：global_sqrt 以接近 human 的误差减少约 60% 耗时；transport_fast 再快约 1.48 倍，但最大相对误差约为 global_sqrt 的 64 倍。若径向目标只要求这批样本中约 `1e-6` 级相对精度，transport_fast 很有吸引力；不能将这个样本结论当成误差保证。

随机 global_sqrt 的最大相对误差出现在案例 983，参考 tau 约 `4.54564e-5`，绝对误差只有 `1.02e-10`。其 P99 相对误差为 `4.58e-8`。因此它的最大相对误差高于 human，并不意味着出现了大光深的严重漏算。三个新增方法的成功结果中均未发现“参考非零却完全算成零”的案例。

为避免把不同成功集合直接比较，六种默认算法共同成功的 **995 个随机案例**如下：

| 算法 | 同一集合耗时（秒） | 最大绝对误差 | 最大相对误差 |
|---|---:|---:|---:|
| human | 0.4032 | 2.33e-7 | 1.55e-6 |
| ref | 1.8585 | 3.51e-7 | 1.59e-6 |
| test_qagp | 0.5801 | 2.31e-7 | 1.55e-6 |
| global_sqrt | **0.2142** | 3.47e-7 | 2.24e-6 |
| velocity_panels | 1.8209 | 1.39e-6 | 1.52e-6 |
| transport_fast | 0.4226 | 1.01e-6 | 9.03e-7 |

transport_fast 默认总耗时与共同成功耗时差别较大，原因是其失败尝试在耗尽求积区间预算前花费了较多时间。失败率和失败成本都应纳入实际选择。

**默认参数与容差对照**

| 算法 | 轨道 atol/rtol | 求积控制 | 最大步长/r |
|---|---|---|---:|
| human、test_qagp | 1e-10 / 1e-10 | atol=rtol=1e-8 | 0.1 |
| ref | 1e-6 / 1e-6 | atol=rtol=1e-6 | 0.01 |
| global_sqrt | 1e-10 / 1e-10 | atol=rtol=1e-8 | 0.5 |
| velocity_panels | 1e-10 / 1e-10 | rtol=1e-6；绝对预算按步长与子区间数分配 | 0.2 |
| transport_fast | 1e-10 / 1e-9 | rtol=1e-7、总绝对预算1e-9，按路径区间分配 | 0.1 |

`matched` 配置把名义轨道 atol/rtol 设为 `1e-10`、求积相对容差设为 `1e-8`。ref 和 transport_fast 的绝对容差/预算设置为 `1e-8`；velocity_panels 保留其 `tolerance*.01*min(1,h/r)/子区间数` 的预算分配。global_sqrt 默认已经使用这些名义轨道及求积容差。不同算法的预算分配与误差估计器不同，这不是严格相同实际误差的比较。

| 样本 | 配置 | 耗时（秒） | 成功/总数 | 最大绝对误差 | 最大相对误差 | P99 相对误差 |
|---|---|---:|---:|---:|---:|---:|
| 径向 | ref_matched | 0.6553 | 900/900 | 4.02e-9 | 3.89e-10 | 1.31e-10 |
| 径向 | velocity_panels_matched | 1.0658 | 900/900 | 6.53e-8 | 5.42e-9 | 3.28e-9 |
| 径向 | transport_fast_matched | 0.1532 | 900/900 | 3.94e-6 | 3.27e-7 | 6.83e-8 |
| 随机非径向 | ref_matched | 1.8830 | 1000/1000 | 2.24e-7 | 1.58e-6 | 1.10e-8 |
| 随机非径向 | velocity_panels_matched | 2.3003 | **1000/1000** | **7.65e-8** | **2.75e-7** | **8.52e-9** |
| 随机非径向 | transport_fast_matched | 2.0366 | 976/1000 | 1.46e-6 | 3.97e-7 | 8.89e-8 |

速度事件分段在收紧求积容差后，非径向精度明显提高。transport_fast 单纯收紧名义容差不能解决所有误差，随机失败从 3 例增加到 24 例。

velocity_panels 的历史驱动器会预计算并复用分位点。额外按这种方式测试，数值结果完全相同：径向积分耗时 **0.9242 秒**，分位点预计算 **0.001285 秒**；随机积分耗时 **1.6961 秒**，分位点预计算 **0.06743 秒**。主表逐条重新准备分布所花的时间分别约为 **0.1270/0.1335 秒**。即使采用历史缓存方式，速度排序也没有改变。

**失败、警告与更严格的审计**

- human：案例 598，`GSL_EROUND`。
- test_qagp：案例 598、996，`GSL_EROUND`。
- transport_fast：案例 334、520、750，均为 E 模，`adaptive quadrature exceeded interval limit`；历史原版逐项复现。
- global_sqrt：案例 574、575、598、599，共 12 次 GSL 初始舍入警告。诊断副本确认，2 次在重试后变为 GSL_SUCCESS，10 次最终仍为 GSL_EROUND，并按原实现的 `error <= 10*tol*(1+|value|)` 条件接受。4 个案例的总光深最大绝对误差为 `2.50e-7`、最大相对误差为 `3.65e-8`。主表的 1000/1000 因此应理解为“按历史策略全部返回”，不能解释为零警告或严格的逐区间 `1e-8` 误差保证。

还对全部 1900 个案例测试了更严格的审计配置。这些只运行一次，不参加主要速度排名：

- global_sqrt：cap、轨道和求积容差都乘 0.25。
- velocity_panels：tolerance=1e-9、cap=0.05r、probes=16。
- transport_fast：历史 audit 设置，轨道容差均1e-11、求积 atol/rtol=1e-11/1e-9、cap=0.025r、禁用尾裁剪、开启磁场表节点分割与极值边界保护。

| 审计配置 | 径向成功数 | 径向最大相对误差 | 非径向成功数 | 非径向最大绝对误差 | 非径向最大相对误差 |
|---|---:|---:|---:|---:|---:|
| global_sqrt_audit | 900/900 | 1.80e-9 | 1000/1000 | 3.08e-6 | 2.24e-6 |
| velocity_panels_audit | 893/900 | 6.84e-10 | 1000/1000 | 7.59e-9 | 5.18e-8 |
| transport_fast_audit | 855/900 | 1.81e-10 | 836/1000 | 6.73e-8 | 2.83e-8 |

误差仍只统计成功结果。global_sqrt 审计的随机警告增加到 40 次，而且最大绝对误差变大；不能假设简单缩小容差必然改善这种跨多个轨道片段和分段线性磁场表的全局求积。velocity_panels 的更严格配置在随机路径上很准，但径向出现 7 次求积不收敛。transport_fast 严格配置成功结果精度较高，失败数量也更多。另测其 `guard_extremum_bounds=true` 的 guarded 模式，默认结果与失败集合均未改变，本批失败不是那个边界保护开关能解决的。

**如何选择**

- 需要同时覆盖径向和一般非径向总光深、重视速度，并接受本批观察到的误差及历史舍入警告策略：优先 **global_sqrt 默认配置**。它在接近 human 的主要误差水平上更快，且所有随机样本均返回。
- 只算径向、允许本批约 `3.5e-7` 的最大相对误差：**transport_fast 默认配置**吞吐最高，比 human 快约 3.66 倍。对一般非径向路径不能直接沿用这一优胜结论。
- 一般非径向精度比速度更重要：**velocity_panels_matched** 是这次重复计时配置中较好的选择，1000/1000 成功，最大绝对/相对误差及 P99 均优于其他 matched 配置，但比 global_sqrt 慢约 10.6 倍。
- 径向高精度：**ref_matched** 的这批结果更好，900/900 成功，最大相对误差 `3.89e-10`；代价是比 global_sqrt 慢约 8.1 倍。

本实验只比较单次自由传播至吸收/逃逸边界的**总光深**。transport_fast 的抽样散射位置反演、提前达到目标 tau 时的优势，以及 global_sqrt 缓存整条轨道的开销，在完整 Monte Carlo 多次散射中可能形成不同排序；本次没有测这些流程。有限随机样本也不能证明事件搜索不存在极窄区间漏检。

**编译、运行和复现**

```sh
make global_sqrt
make velocity_panels
make transport_fast
make benchmark_optical_depth
```

前三条编译对应源码并运行 900 个径向案例，输出到 `output/<算法名>.txt`。最后一条调用 `python3 py/compare_all_optical_depth.py`，执行全部主比较、独立径向参考、细步长参考、重点加密、审计、原程序与历史提交一致性验证。依赖本地三个历史提交、`table/bfield_t10.txt`、`table/bench_od.txt`、g++-16 和 Homebrew GSL，Python 只用标准库。

新程序也可直接接受任意路径文件：

```sh
build/global_sqrt build/optical_compare/random.dat output/global_custom.dat default
build/velocity_panels build/optical_compare/random.dat output/velocity_custom.dat matched
build/transport_fast build/optical_compare/random.dat output/fast_custom.dat guarded
```

输入列：`id beta0 muz alpha az omega_inf r0 pol`，pol=1 为 E，0 为 O。批量结果前九列兼容上一轮：`id tau seconds failure geometry density steps quadrature warnings`；新程序第十列为本次调用中分布准备耗时。不同历史代码的工作量计数口径并不完全相同：velocity_panels 的 quadrature 是节点计算数，transport_fast 是 GK 面板数，global_sqrt 是求积调用数；velocity_panels/transport_fast 失败时未返回内部部分工作量，不能直接将这些计数相加比较算法效率。

完整产物保留在 `build/optical_compare/`：

- `manifest.json`：提交、源码/表格 SHA256、编译参数、配置和随机种子。
- `metrics.json`：逐算法总体、共同成功集合、偏振及表面/磁层子组统计。
- `audit_metrics.json`、`refinement_metrics.json`、`velocity_metrics.json`：附加审计与参考收敛验证。
- `port_validation.json`、`original_validation.json`、`standalone_validation.json`：历史版、原程序及独立入口一致性。
- `radial.dat`、`random.dat`、每轮 `.dat/.log`：全部输入、逐案例输出、失败日志。
- `global_warning_details.log` 与 `global_diagnostics.cpp`：12 次警告的最终状态与接受误差估计。

上一轮 `build/qagp_compare/` 及其报告保留，本报告使用本次统一重测的数据。
