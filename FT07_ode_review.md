**FT07 / resonant.typ / ode 数值实现审查**

审查日期：2026-09-22。依据优先级为 `resonant.typ` > `article/FT07.pdf`。检查了 `ode/` 下全部 13 个 `.cpp` / `.hpp` 文件，并辅助检查了磁场表、生成脚本、Makefile 和已有测地线测试。此次没有修改数值计算源码。

结论：在正常、非退化状态下，核心物理公式基本正确实现了 note 的“单种电子、全部电流由该种粒子承担”的特例。当前程序还不能作为经过可靠性验证的参数扫描和论文统计复现程序：已复现漏掉整个共振层、合法几何状态产生 NaN 后积分器不终止、步端事件漏检等问题。以下明确区分已复现缺陷、模型选择和进一步建议。

**1. 高优先级：仅靠当前 RK 局部误差控制，会漏掉整个共振层。**

位置：[main.cpp:163](/Users/tz_mbp/Desktop/magmc_tz/ode/main.cpp:163)、[dopr5.hpp:98](/Users/tz_mbp/Desktop/magmc_tz/ode/dopr5.hpp:98)、[dopr5.hpp:124](/Users/tz_mbp/Desktop/magmc_tz/ode/dopr5.hpp:124)。

目前没有限制一个步长内共振速度、判别式或粒子分布的变化。当各 RK 采样点都落在散射率很小的位置时，嵌入误差估计也很小，即使采样点之间包含整个高光深共振层。

复现实验：沿径向从 `r=10 km` 积到 `r≈10000 km`，`mu_z=-0.2`、`omega_inf=1 keV`，磁场、质量、半径保持项目现值；把合法温度参数设为 `beta0=-0.2`。使用项目的实际 `PhotonEvolution` 和 `StepperDopr5`，`atol=rtol=1e-6`，初始步长 `0.01 km`，残余光深初值 `-1`。为比较总光深，此实验不执行散射，也不在光深过零处停止。

| 模式 | 原积分器得到的累计光深 | 独立积分参考值 |
|---|---:|---:|
| O | 约 `1.20e-14` | `0.461715239204475` |
| E | 约 `3.42e-14` | `4.66235008107873` |

E 模这条未改变方向的光线，参考光深对应约 `99.1%` 的至少一次散射概率，原结果却相当于没有散射层。这不是 Monte Carlo 抽样涨落。

参考值通过改用粒子速度积分得到，消去了共振 Jacobian，不使用项目的二次方程求根或 RK 光深积分。径向光线的 `mu=k·Bhat` 不变，令 `q=2+p`、`L(r)=sqrt(1-rs/r)`、`Q(r)=q-rs/[2(r-rs)]`，则

\[
\tau=\frac{(p+1)\pi(B_\phi/B_\theta)}{|\bar\beta|}
\int f(\beta)|e'(\beta)|^2\frac{1-\beta\mu}{Q[r(\beta)]L[r(\beta)]}\,d\beta.
\]

`r(beta)` 用单调二分求解共振条件，只计传播区间内的根。20,000 与 40,000 分段的 Simpson 积分结果相符到所列数字。这个公式是本次对 note 的径向特例所作的变量替换，不是直接引用 FT07 的个体速度分母。

对该例限制 `h<=0.01*r` 后能找到共振层；进一步收紧容差可以改善积分精度。但固定相对步长上限本身也不保证解析任意低温分布，不能把这个经验值作为普遍解决办法。

建议：根据共振速度相对于分布宽度的变化限制步长，并显式寻找 `D=0` 和分布支持边界；对 E 模可积奇点采用变量替换或分段积分。note 允许不直接计算精确 `D=0` 点，这没有问题；必须正确积分它的邻域。FT07 §3.3–3.4 的共振速度步长限制值得参考，不必照搬其光深公式。

**2. 高优先级：非物理解没有先过滤，可使实际光深变成 NaN，继而无限循环。**

位置：[main.cpp:60](/Users/tz_mbp/Desktop/magmc_tz/ode/main.cpp:60)、[main.cpp:94](/Users/tz_mbp/Desktop/magmc_tz/ode/main.cpp:94)、[dopr5.hpp:98](/Users/tz_mbp/Desktop/magmc_tz/ode/dopr5.hpp:98)、[dopr5.hpp:144](/Users/tz_mbp/Desktop/magmc_tz/ode/dopr5.hpp:144)。

note 第 241 行要求只保留 `-1<beta<1` 且位于分布支持内的根。代码先计算 `fb.f(beta)`，但即使结果为零，仍然计算 aberration 和含除法的权重；`0*NaN` 不等于零。

已用实际磁场表和默认 `fb(-0.75)` 复现：`r=100 km`、磁赤道、光子沿局部磁场、`x=omega_c/omega=2`。二次方程给出 `beta=-0.6` 与 `beta=1`。后一根必须排除，却产生 `0/0`，使 E/O 两种模式的 `calc_dtaudl` 都返回 NaN。这是在磁赤道发生的缺陷，与下面的磁轴奇点不同。

NaN 误差进入 `success()` 后，步长也变成 NaN；下溢判断对 NaN 返回 false，`do_step()` 的重试循环无法退出。给导数函数加 100 次调用保护后，实际观察到 `h=nan, y_new=nan`。

建议：在任何根相关除法前检查有限性、严格物理速度范围和分布支持；对零分布值直接跳过。积分器也应拒绝非有限导数、状态、误差和步长，输出光子状态后终止该条轨迹。重建的 `mu` 在舍入误差内越界时可以钳制到 `[-1,1]`，明显越界则应报错。

**3. 已复现：恰好落在步端的事件会漏掉。**

位置：[dopr5.hpp:183](/Users/tz_mbp/Desktop/magmc_tz/ode/dopr5.hpp:183)。

条件仅接受严格异号。若新端点的事件值恰为零，这一步不触发事件，随后记录的符号又是零，下一步也不触发。

最小复现：`y'=1`，事件函数为 `x-1`，`init(0,1,{0})`。首步到 `x=1` 返回 `event=-1`；下一步到 `x=11` 仍返回 `event=-1`。逃逸、吸收和光深阈值都有相同问题。

应明确处理“从非零到达零”的端点，同时区分从事件面出发的方向。不能直接把所有起点零都当吸收：程序最初就在星面向外发射。

**4. 已复现：磁轴处磁场和局部基底出现 0/0。**

位置：[bfield.hpp:42](/Users/tz_mbp/Desktop/magmc_tz/ode/bfield.hpp:42)、[main.cpp:52](/Users/tz_mbp/Desktop/magmc_tz/ode/main.cpp:52)、[main.cpp:67](/Users/tz_mbp/Desktop/magmc_tz/ode/main.cpp:67)、[main.cpp:84](/Users/tz_mbp/Desktop/magmc_tz/ode/main.cpp:84)。

`calc_B(10,+1)` 实际返回 `(1e14,NaN,NaN)`，`calc_B(10,-1)` 返回 `(-1e14,NaN,NaN)`；正确极限分别是 `(±1e14,0,0)`。此外 `rho=0` 时 `theta_hat/phi_hat` 的构造也除零。

应实现极轴解析极限，并用稳定的局部笛卡尔场方向处理坐标退化。电流系数直接计算 `Bphi/Btheta=A*f^(1/p)`，不要在两分量同时趋零时作除法。只修 `calc_B` 还不能修复整条计算链。对恰好径向的出射光子，`to_unit(cross(rhat,kout))` 同样需要任意垂直平面的后备选择。

**5. 已复现：合法的低温 Boltzmann 参数发生共同欠流。**

位置：[distribution.hpp:9](/Users/tz_mbp/Desktop/magmc_tz/ode/distribution.hpp:9)、[distribution.hpp:22](/Users/tz_mbp/Desktop/magmc_tz/ode/distribution.hpp:22)。

`beta0=-0.05` 满足 note 的参数条件，但此时 `a≈798.50`，`K1(a)` 和 `exp(-a)` 均欠流为零，均值及概率密度变成 NaN。默认 `beta0=-0.75` 不受此问题影响。

建议使用真正的缩放 Bessel 函数 `K1e(a)=exp(a)*K1(a)` 或其渐近式/对数形式：

\[
f(\beta)=\frac{\gamma^3e^{-a(\gamma-1)}}{K_{1e}(a)},\qquad
\bar\beta=\frac{\operatorname{sgn}(\beta_0)}{aK_{1e}(a)}.
\]

不能先算已欠流的 `K1` 再乘 `exp(a)`。还应验证 `0<abs(beta0)<1`，并避免小 `beta0` 下直接计算 `gamma0-1` 的消减误差。`f(0)` 的当前处理和 note 的严格支持定义也不同；单点本身不改变连续分布积分，但实现支持边界时宜保持一致。

**6. 确定的保护缺失：零总权重仍会执行散射。**

位置：[main.cpp:102](/Users/tz_mbp/Desktop/magmc_tz/ode/main.cpp:102)。

有两个实根不等于有可散射粒子。若两根都不在分布支持内，总权重为零，代码仍会选第二根并更新光子。总权重为 NaN 时也会走到类似分支。当前 4,000 条默认参数试运行没有触发这一情况，因此不声称默认运行已产生这种伪事件。

应要求总权重有限且严格大于零；若事件插值给出零光深处的散射点，应回退并重新定位。建议光深与选枝共用同一个“有效根及权重”函数，避免两份实现以后不一致。O 模单根权重还可稳定地化简为 `0.5*f(beta)*(1-beta^2)*abs(mu-beta)`。

**7. 磁场插值误差不能由 ODE 容差控制。**

位置：[bfield.hpp:39](/Users/tz_mbp/Desktop/magmc_tz/ode/bfield.hpp:39)。

磁场表本身通过了边界、BVP 和扭转角核对。但是当前分别线性插值 `f` 和 `f′`，两者不再严格满足导数关系，也不严格保留解析场的无散性质。当前 `Delta_mu=0.001`，最后一个网格区间靠近极点时，`Btheta` 相对正确极区展开的误差趋向约 `-4.17e-4`，`Bphi/Btheta` 约 `-4.71e-4`。

这不表示最终谱一定有同样大小的误差，但 `rtol=1e-6` 不能代表全模型达到该精度。建议用同一 Hermite 插值多项式生成 `f` 和它的导数，并处理极区渐近式；至少比较磁场表加密前后的可观测量。

**8. 与论文的模型差异和目前实现范围。**

- [main.cpp:184](/Users/tz_mbp/Desktop/magmc_tz/ode/main.cpp:184) 将所有返星光子吸收。FT07 §3.3，期刊第 622 页，采用 E 模全吸收、O 模半数弹性表面散射、半数吸收。note 没有指定表面反照率，因此这应视为需要明确的模型选择；若目标是相同设置的 FT07 复现，应补齐该边界处理。参见[论文](https://arxiv.org/abs/astro-ph/0608281)。
- 当前只有一类电子、`epsilon=1`。这正确实现 note 的特例，但尚未实现多种粒子分别设定电流份额、质量和分布后求和/选种类的通用形式。
- [main.cpp:148](/Users/tz_mbp/Desktop/magmc_tz/ode/main.cpp:148) 只演化一个光子；还没有集合循环、输出能量/观测方向统计、谱和脉冲轮廓。逃逸时打印的 `psi,alpha` 是该传播段平面内变量，不能直接当成观测者相对磁轴的方向。
- [init.hpp:15](/Users/tz_mbp/Desktop/magmc_tz/ode/init.hpp:15) 的径向、均匀表面、单能种子是有效的特定模型设置，本身不是错误。`omega_inf=1` 配合当前常数实际表示远处能量 `1 keV`，表面局部能量约为 `1.306 keV`。若希望表面为 `1 keV`，初始化应乘表面 lapse。
- `r=1000*R_star` 是有限外边界。应通过增大外边界验证残余散射和观测方向的收敛；还应从最终平面基底重建方向，必要时继续真空测地线到渐近方向。
- 引入 Schwarzschild 传播并修正 FT07 的密度归一化后，即使代码完全正确，也不保证逐点重现论文原图。后续比较应说明这些模型差别。GR 加固定半径也使直接假定原论文的频率缩放/单一响应核需要重新验证。

**已确认与 note 一致的关键部分。**

- 三个 Schwarzschild 光线方程、局部频率红移和光子方向余弦的坐标变换正确。`r,l` 同用 km 时，光深公式的 `1/r` 也给出每 km 的量，不需要另外补一个长度单位换算。
- 速度密度中的 `gamma^3` Jacobian 正确。分布的归一化和均值公式与 note、[NIST DLMF 10.32.9](https://dlmf.nist.gov/10.32.E9) 一致。默认 `beta0=-0.75` 的平均速度是约 `-0.48802144336328246`，不是 `-0.75`。
- **光深前因子 `1/abs(fb.b_bar())` 正确，不应改回 FT07 印刷式里的逐根 `1/abs(beta)`。** 这里采用的是 note 的数目加权分布和电流密度定义。
- 非退化物理根上的共振方程、逐根 O/E overlap、选枝权重、出射 aberration、能量比、偏振概率和轨道平面更新均一致。
- 南北半球延拓正确：`Br` 变号，`Btheta/Bphi` 保持相应对称性；正扭转约定一致。
- `sample_mup` 的实现虽然不像直接逆 CDF，但严格给出 `3*(1+mu_prime^2)/8`。它是权重 `3/4` 的均匀分布与权重 `1/4` 的对称三次幂分布的混合。
- DOP5 的阶段、嵌入误差、dense output 系数没有发现抄写错误；稳定二次求根的缩放与 q 公式方向正确。
- `B_to_omega` 实际把磁场换成 `hbar*omega_c` 的 keV 数值。只要全部 `omega` 都采用相同能量表示，当前比例运算自洽；建议把命名和注释改得更明确。

**较低优先级的健壮性与维护建议。**

- [ran.hpp:35](/Users/tz_mbp/Desktop/magmc_tz/ode/ran.hpp:35)：整数转换为 double 后可得到端点 `U=0` 或 `U=1`，而指数自由程需要开区间均匀数。`log(0)` 为负无穷，`log(1)` 产生零阈值。发生概率极低，但宜明确排除；圆周拒绝采样也应排除 `r2=0`。
- [roots.hpp:17](/Users/tz_mbp/Desktop/magmc_tz/ode/roots.hpp:17)：未缩放的平方/乘积会欠流。实际 `zriddr(1e-200*(x-0.2),0,1,1e-12)` 返回初始哨兵 `-9.99e99`，而非 `0.2`。当前事件通常不处于这个量级，但基础求根器应修复退化分支或采用缩放。
- [solve_quadratic.hpp:1](/Users/tz_mbp/Desktop/magmc_tz/ode/solve_quadratic.hpp:1) 缺少 include guard；重复包含会报重定义。共振根刚合并时应以专门的共振条件和残差验证辅助一般二次求根器，而不是把其所有失败都等同于物理无共振。
- `BField` 读取文件后未验证表头和每行提取成功；`UniformHunt` 未验证节点数、区间和有限输入。有效现有表下未发现索引错误，但坏表可能成为静默物理错误或非法内存访问。
- [test/test_null_geodesic.cpp:28](/Users/tz_mbp/Desktop/magmc_tz/test/test_null_geodesic.cpp:28) 仍调用旧版七参数构造函数；即使传 `-Iode` 也不能编译，需更新为构造后调用 `init()`。现有测试不能直接作为本版积分器验证。

**本次执行的验证及其限度。**

| 检查 | 结果 |
|---|---|
| GCC 16，C++23，`-O2 -Wall -Wextra -Wpedantic` 编译主程序 | 成功；只有未使用参数警告 |
| 原始种子 1234、O 模、默认物理参数 | 在 `r≈121.169981 km` 散射一次，随后逃逸 |
| 默认物理参数，种子 0–1999，O/E 各 2000 条 | 3958 条逃逸，42 条吸收，3399 次散射；未触发测试的非有限值/非法散射保护 |
| 低温参数的径向光深，与独立速度积分比较 | 复现严重漏共振；见第 1 项 |
| 默认 `beta0=-0.75` 的若干径向总光深收敛 | `1e-6` 容差仍可有数 `1e-4` 的总光深误差；收紧后改善，不应把局部容差等同于全局误差 |
| 无散射测地线，初始 `r=100 km`、`alpha=0.5,pi/2,2.5` | 冲量参数 `r*sin(alpha)/sqrt(1-rs/r)` 最大相对漂移：容差 `1e-6` 时约 `9.2e-8`，`1e-9` 时约 `2.0e-10` |
| `beta0=±0.5,±0.75` 分布数值积分 | 归一化误差小于约 `2.1e-14`，均值误差小于约 `3.3e-13` |
| 百万次 `sample_mup` | 二阶矩 `0.400369`，理论 `0.4`；四阶矩 `0.257471`，理论 `9/35` |
| 磁场表 | `C=A^2*p*(p+1)`、边界值一致；离散 BVP 最大相对残差约 `1.32e-7`；重积分扭转角 `1.000000000175 rad` |
| 极轴、平行场方向、事件步端、NaN 重试、低温欠流 | 均有确定复现，见上文 |

4,000 条正常结束仅证明这组试运行没有触发所加的保护，不等于证明光深、光谱或稀有路径均无偏。已有明确反例足以说明还需要数值修复和收敛验证。

临时复现源文件保存在 `/tmp/magmc_review/`，包括 `radial.cpp`、`radial_focus_residual.cpp`、`ensemble.cpp`、`geodesic.cpp`；其他边界测试位于 `/tmp/magmc_aligned_probe.cpp`、`/tmp/magmc_integrator_test.cpp`、`/tmp/magmc_nan_test.cpp`、`/tmp/ft07_field_sampling_review.cpp`。这些是审查用临时文件，系统清理临时目录后可能消失。例如从项目根目录运行：

```sh
g++-16 -std=c++23 -O2 /tmp/magmc_review/radial_focus_residual.cpp -o /tmp/magmc_review/radial_focus_residual
/tmp/magmc_review/radial_focus_residual
```

这些 harness 通过重命名并包含主程序来调用实际实现；编译时可能出现被重命名的原 `main` 缺少显式返回值警告。该函数不被测试调用，测试入口为另一个 `main`。

**逐文件覆盖记录。**

| 文件 | 审查结论 |
|---|---|
| `main.cpp` | 正常根物理公式正确；根筛选、散射保护、边界和统计驱动有上述问题/缺项 |
| `bfield.hpp` | 公式、南北对称性正确；极轴、插值和读表验证需改进 |
| `constants.hpp` | 当前单位运算自洽；建议明确 keV 表示的是能量 |
| `distribution.hpp` | 正常温度公式正确；低温欠流及参数/支持验证需改进 |
| `dopr5.hpp` | 核心系数正确；端点事件、非有限值和共振步长处理有问题 |
| `doubles.hpp` | 向量代数正确；零向量归一化须由函数或调用处处理 |
| `hunt.hpp` | 当前有效表的索引正确；需输入验证 |
| `init.hpp` | 特定种子模型有效；应明确能量位置及初始偏振设置 |
| `photon.hpp` | 状态内容与 note 一致，未发现独立缺陷 |
| `ran.hpp` | 球面、圆周、垂直方向抽样正常；浮点端点和退化输入需处理 |
| `roots.hpp` | 常规括根用途正常；极小函数尺度可错误返回哨兵 |
| `sample_mup.hpp` | 概率分布正确，有解析与抽样验证 |
| `solve_quadratic.hpp` | 稳定求根框架正确；物理根过滤须由调用层完成，并补 include guard |

建议修复顺序：先避免漏共振及 NaN 不终止，再处理事件/坐标/分布边界，之后完成磁场表与误差容差的收敛测试，最后扩展为记录能量、渐近观测方向、偏振和散射次数的光子集合计算。
