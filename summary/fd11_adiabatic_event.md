# 用绝热条件触发偏振积分的测地线事件

本次依据用户提供的 `doc/adiabatic.md` 修改 `src/dev/fd11.cpp`。启动条件改为

\[
\left|\frac{d\chi_B}{dl}\right|>\epsilon_{\rm ad}\kappa,
\qquad \epsilon_{\rm ad}=10^{-3}.
\]

文档用复振幅表达该条件，实际程序仍按此前约定演化三个实数 Stokes 参数。
`--couple` 现在指定这里的 $\epsilon_{\rm ad}$；默认数值仍是 `1e-3`，不再表示旧的 $1/(\kappa r)$ 阈值。
冻结继续使用原有 FD11 (35) 判据。

## 1. 从随磁场转动的基底得到启动事件

局域本征模振幅满足用户文档给出的方程

\[
\frac{d}{dl}\begin{pmatrix}a_O\\a_E\end{pmatrix}
=\begin{pmatrix}i\kappa/2&\chi'_B\\-\chi'_B&-i\kappa/2\end{pmatrix}
\begin{pmatrix}a_O\\a_E\end{pmatrix}.
\]

两个对角传播本征值的差是 $i\kappa$，非对角耦合项的大小是 $|\chi'_B|$。
所以比较的无量纲量为 $|\chi'_B|/\kappa$，而不是 $|2\chi'_B|/\kappa$。
当这个比值小于阈值时仍记录 O/E 标签；达到过渡面时初始化 Stokes 并开始数值演化。

为避免显式除以可能很小的 $\kappa$，定义事件函数

\[
G(l)=|\chi'_B(l)|-\epsilon_{\rm ad}\kappa(l).
\]

两项的单位都是 `km^-1`。绝热区中 $G<0$，非绝热区中 $G>0$，共享 stepper 定位 $G=0$。
事件根本身作为过渡位置；初始位置恰好位于过渡面时也按开始积分处理。

## 2. 如何计算角度导数

事件在 `stepper.init()` 时也会求值，这时轨道稠密输出尚未准备好。
因此 `chi_prime(state)` 直接对给定几何状态求局部方向导数，不依赖 `dense_out` 或旧的轨道步长。

几何状态为 $y=(r,\psi,\alpha)$，沿局域路程的切向量由现有测地线方程给出：

\[
v=\frac{dy}{dl}=\left(
\mathcal N\cos\alpha,
\frac{\sin\alpha}{r},
-\frac{\sin\alpha}{r\mathcal N}\left(1-\frac{3r_s}{2r}\right)
\right),\qquad \mathcal N=\sqrt{1-r_s/r}.
\]

取 $\delta=10^{-5}r$，构造局部差分采样点 $y_-=y-\delta v$、$y_+=y+\delta v$。
对一个光滑状态函数 $f$，Taylor 展开为

\[
f(y\pm\delta v)=f(y)\pm\delta\nabla f\cdot v
+\frac{\delta^2}{2}v^T(\nabla\nabla f)v+O(\delta^3).
\]

两式相减并除以 $2\delta$ 后，偶次项抵消：

\[
\frac{f(y_+)-f(y_-)}{2\delta}=\nabla f\cdot v+O(\delta^2)
=\frac{df}{dl}+O(\delta^2).
\]

现有 `pol_coeff` 返回 $c=\cos(2\chi_B)$、$s=\sin(2\chi_B)$。
记 $\Theta=2\chi_B$，则两个采样点之间有

\[
c_-s_+-s_-c_+=\sin(\Theta_+-\Theta_-),
\qquad c_-c_++s_-s_+=\cos(\Theta_+-\Theta_-).
\]

因此先计算局部角差

\[
\Delta\Theta=\operatorname{atan2}
(c_-s_+-s_-c_+,\ c_-c_++s_-s_+).
\]

它避免直接相减带分支跳跃的 `atan2` 角度。采样点相距 $2\delta$，同时 $\Theta'=2\chi'_B$，所以

\[
\boxed{\chi'_B\simeq\frac{\Delta\Theta}{4\delta}}.
\]

分母中的两个 2 分别来自采样间距和双角关系，不能遗漏其中任何一个。

对所有路径统一将 $|\Delta\Theta|\le16\epsilon_{\rm machine}$ 的角差视为数值上不可分辨的零。
否则，极弱场中舍入噪声本身就可能超过 $10^{-3}\kappa$。这不是径向光子的特殊分支。
另外，事件回调在当前轨道步右端显式使用 `stepper.y_new`，使端点变号检查与求根时的端点值一致；
稠密输出在端点的末位舍入差异不再破坏这个事件的括根条件。

## 3. 事件和状态流转

事件注册顺序为：`AbsorbEvent=0`、`EscapeEvent=1`、`PolarizationEvent=2`。
普通路径仍通过 `stepper.detect_event()` 选择最先到达的事件，并截断轨道步。

`evolve_geodesic()` 的顺序现在是：

1. `step_geodesic()` 推进轨道并取得事件编号。
2. `advance_polarization()` 仅在状态已经为 `Integrating` 时推进偏振；`Mode` 状态直接返回。
3. 对截断后的当前区间计算光学深度、处理可能更早发生的散射。
4. 若未被更早的散射中断，在函数循环末尾处理事件编号；第三个事件调用 `start_polarization()`。
5. 在事件点设置 `pol_start`，用局域本征模初始化 Stokes，切换到 `Integrating`，然后继续下一段轨道。

因此，若散射先于原本的启动事件发生，不会提前把该光子切换为积分状态。
事件不重新抽样光学深度，已累计的光学深度继续保留。

启动后的本征模初值仍为

\[
(Q,U,V)_O=(\cos2\chi_B,\sin2\chi_B,0),\qquad
(Q,U,V)_E=-(\cos2\chi_B,\sin2\chi_B,0).
\]

对于发射或散射后已经满足条件的光子，`step_geodesic()` 根据初始化时记录的事件符号，
提供当前位置的零长度 `PolarizationEvent`；同样由循环末尾的事件分支启动。
这避免只检测变号时漏掉一开始就在非绝热区内的光子。

进入 `Integrating` 后，启动事件回调返回固定负值，并同步重置该事件的缓存符号，避免禁用时制造虚假过零。
冻结后也不再触发启动事件；下一次散射重新初始化为 O/E 模式时，通过 `stepper.init()` 重新计算事件符号。
关闭偏振演化时，启动事件一直返回负值。

原来的式 (34) 判据、四段启动扫描及 `advance_polarization()` 内的独立求根已经删除。
光学深度中用于共振边界和粒子分布节点的扫描仍保留。启动事件沿用共享 stepper 的步端变号检测，
不额外保证捕捉同一轨道步内先越界又返回的任意窄区间；最大轨道步长仍是 `0.1*r`。

## 4. 验证与误差说明

`make fd11-test` 通过。测试覆盖新事件注册、角度导数的两倍关系、延迟启动、初始即非绝热、
单次触发、散射后重新启用、关闭偏振、常方向纯模式，以及原有代数、散射率和分部积分检查。

独立 Jones 参考解测试得到：

|初始传播角|启动半径 km|默认 `pol-tol=1e-5` 最大 Stokes 差|`pol-tol=1e-8` 最大 Stokes 差|
|---:|---:|---:|---:|
|0.4|398.976210|6.163e-7|3.226e-7|
|1.0|180.000000（初始点）|1.784e-5|1.408e-6|
|1.8|175.312822|6.689e-4|5.500e-7|

新判据确实改变了积分起点。最后一条近磁场平行轨道在默认容差下超过了旧测试使用的 `2e-4` 全局误差上限，
收紧容差后明显收敛。当前回归分别要求默认轨道误差小于 `1e-3`、收紧容差后的独立 Jones 误差小于 `2e-5`，
没有通过改变生产代码的默认 ODE 容差来掩盖这一差异。局部容差和最终全局误差不是同一个量。
其中 `alpha=1.8` 初始向内传播，旧启动半径约 `176.928 km`，新启动半径约 `175.313 km`，
即沿该光路稍晚启动；不能仅由误差增大推断它更早启动。

此外完成八组端到端模拟，每组 1000 光子，种子 4193：

|初始模式|能量 keV|偏振演化|逃逸|吸收|启动次数|逃逸时仍为 Mode|
|---|---:|---|---:|---:|---:|---:|
|E|0.1|开启|1000|0|739|532|
|E|1|开启|992|8|476|538|
|E|10|开启|963|37|427|537|
|O|0.1|开启|996|4|477|686|
|O|1|开启|992|8|309|690|
|O|10|开启|959|41|281|678|
|E|1|关闭|992|8|0|992|
|O|1|关闭|992|8|0|992|

启动次数按每段自由传播累计，发生散射后可以重新启动。所有输出均为有限的五列数据，
八组中的最大 `|Q²+U²+V²-1|` 为 `2.52e-13`。
当前程序仍采用随机表面位置、初始径向发射的原驱动；投影磁场方向不转动时，
新条件自然不会触发数值积分，所以部分光子以 Mode 状态逃逸是预期行为。

回归日志为 `output/fd11_adiabatic_test.log`；模拟记录为 `output/fd11_adiabatic_validation.json`，
各组光子数据及日志位于同目录的 `fd11_adiabatic_*.txt` 和 `fd11_adiabatic_*.log`。
编译产物仍位于 `build/`，`git diff --check` 通过。

## 5. 复现

```sh
make fd11-test
build/fd11 --photons 1000 --seed 4193 --energy 0.1 --output output/fd11_adiabatic_E_0.1_pol.txt
```

需要更严格的轨道偏振精度时，可保留新的启动判据并设置 `--pol-tol 1e-8`；
这与改变 `--couple` 所定义的物理初始化近似是不同的调整。
