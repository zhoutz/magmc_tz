# 任意方向光线的共振积分

本次移除了原先依赖向外径向单调性的 `RadialResonance`、`RadialOpticalDepth`、径向速度积分和径向温度限制器，包括对应 benchmark、二进制及过时结果。`py/bench_od.py` 和 `table/bench_od.txt` 是用户指定的参考，保留不动。新的输运算法不反解共振半径，不固定光子–磁场夹角，也不要求半径单调。

## 通用状态和物理量

`transport::integrate(field, distribution, photon, options)` 接受 `Photon` 中的轨道基矢 `n,e1,e2`、初始 `r,psi,alpha`、能量和偏振。它沿局域静止观测者测量的路径长度 l 演化

\[
dr/dl=L\cos\alpha,\quad d\psi/dl=\sin\alpha/r,\quad
d\alpha/dl=-\frac{\sin\alpha}{rL}\left(1-\frac{3r_s}{2r}\right).
\]

每个采样点重新构造磁场方向和光子方向，求 `x=ωc/ω`、`μ=k̂·B̂`。红移、角度变化、磁场角结构均包含在内。南北半球、正负速度分布和两条可贡献的速度根使用同一条计算路径。

`make_photon(r, muz, alpha, azimuth, energy, pol)` 只是初始化工具；azimuth 是局域切平面内相对 θ̂ 的方位角。散射后已有轨道基矢时，可以直接构造 `Photon`，不必使用此工具。

对于 `D=x²+μ²−1>0`，保留速度分布支撑内的所有根。用 Vieta 关系稳定求根，并将 note 中的因子改写为

\[
|\mu-\beta|=\frac{\sqrt D(1-\beta\mu)}x,\qquad
P_E=1/2,\quad P_O=D/(2x^2).
\]

这是数值等价变形，不改变速度分布的平均速度归一化。`expm1(2 logx)` 用于稳定计算 `x²−1`。3000 个任意角度点将新公式与原 `PhotonEvolution::calc_dtaudl` 交叉比较。

## 共同的边界和奇点处理

三种修复方案共用一个通用的轨道步内积分器，比较的是步长控制策略：

1. DOPRI5 只演化三维几何状态，并提供连续插值。几何容差取 `min(1e-10, 0.01*tau_tolerance)`。
2. 在每个已接受的轨道步内定位 `D=0`、β=0，以及多个预设速度对应的共振面
   `Hβ(l)=log x(l)−log(1−βμ(l))+½ log(1−β²)=0`。
3. 速度标记取 β₀ 的若干倍和 relativistic tail 的若干值；它们仅用于分段，完整分布的尾部不丢弃。FT07 对照另有明确的 99.8% 截断。
4. 不只检查步的两个端点。默认使用 8 个子区间探测表面变化；检测到趋势反转时，用黄金分割搜索内部极值，再将极值加入括根节点。因此，即使起点、终点以及普通采样点都在层外，也能发现内部的窄共振区间。局部根求解容差为步内参数的 `3e-15`。
5. 分段后，在每一段 `[l_a,l_b]` 内作通用的路径长度变换
   `l=l_a+(l_b−l_a)sin²(πz/2)`，再使用开节点 8/16 点 Gauss–Legendre 自适应积分。这同时处理两端可能出现的 E 模平方根奇点。变换的是 **l，不是 r**，所以转向、非径向传播和分支合并不破坏它。
6. 距离已定位的 D=0 边界极近时，用单侧 Taylor 系数锚定 D=0，避免浮点消减放大为虚假的奇点误差。该局部处理仍使用真实轨道上的 x、μ 及其变化。

在轨道步内定位 stellar surface 和 escape surface；另外检查径向转向处的最小半径，避免一步的两端都在星表外，却在中间穿入星体。返回结果区分逃逸、撞击星表、达到指定累计光深和达到路径长度上限。

这种发现边界的方法是数值策略，不是对任意未解析高频磁场的数学证明。当前自相似场上的系统角度测试、加密检查及独立窄区间测试见报告。精确退化的高阶切触不应当由普通简单根误差估计来认证；本次不声称有限测试穷尽连续参数空间。

## FT07 控制器

代码：`benchmark/ft07.cpp`，控制器在 `src/transport.hpp`。

实现论文 §3.3–3.4 的步骤，分布边界取 CDF 的 0.001、0.999 分位数，中心取中位数。层外上限 r/10；接近边缘时将其解析到相应半宽的 1/100；层内限制

\[
|\Delta\beta|\le0.01(1-\beta^2)|\beta|.
\]

与上次实现不同，现在使用完整方向导数。由共振条件得

\[
\frac{d\beta}{dl}=\frac{1-\beta^2}{\beta-\mu}
\left[(1-\beta\mu)\frac{d\ln x}{dl}+\beta\frac{d\mu}{dl}\right].
\]

式 (38) 的限制也保留 `dμ/dl`，即 `Δl≤0.01|μ−β|/|dμ/dl−dβ/dl|`。方向导数由沿完整 geodesic 的对称差分计算，差分长度为 `1e-5 r`。

该对照保留论文中央 99.8% 粒子截断，并非完整分布的高精度算法。`D≤1e-6` 的根合并邻域交给共同的边界/奇点积分器，否则式 (38) 会使总光深积分无限逼近边界而不穿过。论文原本在 Monte Carlo 散射后改变光子方向，本任务需要未散射轨道的整个光深，因此明确采用这个数值补充。

因此本次名称的准确含义是 **FT07 控制器 + 通用边界/奇点处理 + GR 几何**，不是原论文完整程序的逐行复刻，也不能将其时间解释为 FT07 原始代码耗时。

## 新策略一：共振变量变化限制器

代码：`benchmark/phase_cap.cpp`。

冷共振层在最不利的角度处，log(x) 尺度可小到 O(β₀²)。同时限制 x 和 μ 的变化：

\[
v=|d\ln x/dl|+\frac{|d\mu/dl|}{\max(0.01,|\beta_0|)},\qquad
\Delta l\le\min\left[r/10,\frac{\eta\beta_0^2+d/2}{v}\right].
\]

其中 d 是当前位置到热分布核心对应共振面范围的相位距离；处在范围内部时 d=0，远离核心时逐渐放宽限制。核心以 `|β|≤min(0.999,4|β₀|)` 内的速度标记估计；这只影响步长，**不截掉区间外粒子的贡献**。零导数附近仍受几何上限约束；共同的步内边界定位负责检查整个被接受步，而不把局部导数近似当作不会漏层的充分证明。

默认 η=0.1，另外测量 0.03 和 0.3。它不需要知道未来共振半径，并且 `dμ/dl` 对任意方位角均有效。

## 新策略二：轨道步内边界定位

代码：`benchmark/event_guard.cpp`。当前 main 使用此方法。

不按热宽度限制每一段空间步长；让几何 DOPRI5 选择步长，并设置 `Δl≤min(0.4,4η)r`，默认 η=0.1。由上述共同积分器定位该步内所有发现的共振面，分段后累计光深。这避免为了穿过一个窄层而在整段轨道上持续取热宽度级别的小步。

这个方法的效率优势来自把几何误差控制与共振积分分开。它没有径向快路径；基准输入虽为径向，执行的仍是同一套任意方向轨道代码。

## 累计光深与散射位置

`options.path_limit` 可指定 proper path length 截止位置；`options.tau_target=-log(U)` 可定位第一次达到抽样光深的位置。发生散射时，只在包含目标的面板内部对累计光深括根，并返回 l 和 `(r,psi,alpha)`。这是自由传播段的定位器；它不代替散射后方向、能量和偏振的采样。

```cpp
Boltzmann charges(-0.3);
auto photon = transport::make_photon(30, -0.4, 2.1, 0.7, 1.0, Polarization::E);
transport::Options options;
options.method = transport::Method::event_guard;
options.tau_target = -std::log(0.5);
auto result = transport::integrate(field, charges, photon, options);
```

这个例子是南半球、非径向且初始向内的光子，使用与 main 和径向 benchmark 相同的求解器。

## 文件与复现

- `src/resonance.hpp`：完整几何下的共振参数、两条速度根与稳定 opacity。
- `src/resonance_panels.hpp`：步内共振面定位、极值探测与奇点积分。
- `src/transport.hpp`：控制器、几何演化、累计光深与事件。
- `benchmark/baseline.cpp`：原四维 DOPRI5 对照。
- `benchmark/{ft07,phase_cap,event_guard}.cpp`：三种通用修复方案。
- `benchmark/angular_validation.cpp`：非径向和正负速度分布验证。
- `benchmark/transport_validation.cpp`：独立窄层测试、opacity 公式、散射事件、轨道守恒量及 FT07 截断分解。
- `output/run_benchmarks.py`：完整构建、计时和验证脚本；其他辅助文件、结果和说明也都位于 output/。

从仓库根目录运行：`python3 output/run_benchmarks.py --repeats 5`。
