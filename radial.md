# 讲解光深积分参考实现

> 我没看懂 "test/resonance\_reference.hpp" 这个文件中，光深积分的参考实现是怎么推导出来的，请帮我用公式讲解一下。

<details><summary>3 previous messages</summary>

> 我先读一下参考实现及其相关定义，再按代码中的变量把光深积分的推导写成公式。
>
> <details><summary>Explored 6 files, a search, a list, ran 2 commands</summary>
>
> - Ran `pwd; rg --files -g 'AGENTS.md' -g '*resonance*' -g '*optical*'`
>
> <details><summary>Explored 2 files</summary>
>
> - Read `./ode/resonance_transport.hpp`
> - Read `./test/resonance_reference.hpp`
>
> </details>
>
> - Ran `rg --files -g AGENTS.md /Users/tz_mbp/Desktop /Users/tz_mbp/.codex 2>/dev/null`
>
> <details><summary>Explored 6 files, a search, a list</summary>
>
> - Listed files
> - Read `./ode/distribution.hpp`
> - Read `./ode/transport_physics.hpp`
> - Read `./ode/bfield.hpp`
> - Searched for `tau|delta|sigma|Boltz|bessel|K_1|gamma|光深|截面|分布|电流` in `resonant.typ`
> - Read `./test/resonance_reference.hpp`
> - Read `./ode/distribution.hpp`
> - Read `./resonant.typ`
>
> </details>
>
> </details>
>
> 关键是交换积分顺序：普通实现先在每个位置求共振速度，参考实现则固定粒子速度，先积掉共振位置。这样既能消去速度根合并时的奇异分母，也能推导出代码里的 \(1/(LQ)\)；随后再换成变量 \(u\)，把低温时很窄的分布展开。
</details>

这个参考实现的核心是：**把“沿光线路径积分”改写成“对粒子速度分布积分”**。先用共振条件消去半径 \(r\)，再把速度 \(\beta\) 换成适合低温分布的变量 \(u\)。

下面推导 [resonance_reference.hpp 中第 79–85 行的两条公式](./test/resonance_reference.hpp:79)。

**先写出还没有消去 \(\delta\) 函数的光深积分。**

按照项目 [resonant.typ 中的定义](./resonant.typ:179)，对代码考虑的单个粒子族，令
\[
C=(p+1)\pi\frac{B_\phi}{B_\theta},
\qquad
\mathcal P(\beta)=|e'|^2,
\]
则
\[
\boxed{
\tau
=
C\int_{\rm path}\frac{dl}{r}
\int_{-1}^{1}d\beta\,
\frac{f(\beta)}{|\bar\beta|}
(1-\beta\mu)\mathcal P(\beta)\,
\omega_D\,\delta(\omega-\omega_D)
}
\tag{1}
\]

这里：

- \(\beta=v/c\)，\(\gamma=(1-\beta^2)^{-1/2}\)；
- \(f(\beta)\) 是相对于 \(d\beta\) 归一化的速度分布，\(\bar\beta=\int\beta f(\beta)d\beta\)；
- \(\mu=\hat{\boldsymbol k}\cdot\hat{\boldsymbol B}\)，是光子方向与磁场方向夹角的余弦；
- \(\omega_D=\omega_c/[\gamma(1-\beta\mu)]\)，共振要求 \(\omega=\omega_D\)。

式中的 \(1/|\bar\beta|\) 来自**用给定电流反推粒子数密度**：相同电流下，平均速度越小，所需粒子越多。

**径向光线使半径依赖变得特别简单。**

这个参考实现只处理 \(\alpha=0\) 或 \(\pi\) 的径向传播。沿途纬度不变，而项目磁场的各分量具有相同的径向因子：
\[
B(r)=B(R_*)\left(\frac{R_*}{r}\right)^q,
\qquad q=2+p.
\]
因此磁场方向、\(B_\phi/B_\theta\) 和 \(\mu\) 都沿路径保持不变。向外、向内分别有
\[
\mu=\pm\frac{B_r}{B}.
\]

定义引力红移因子
\[
L(r)=\sqrt{1-\frac{r_s}{r}}.
\]
局域光子频率和径向局域路程满足
\[
\omega(r)=\frac{\omega_\infty}{L(r)},
\qquad
dl=\frac{|dr|}{L(r)}.
\]

再定义
\[
x(r)=\frac{\omega_c(r)}{\omega(r)}
=
\frac{\omega_c(R_*)}{\omega_\infty}
\left(\frac{R_*}{r}\right)^qL(r),
\qquad
g(\beta)=\gamma(1-\beta\mu).
\]
于是共振条件就是
\[
\boxed{x(r)=g(\beta).}
\tag{2}
\]

普通实现固定 \(r\)，求出可能的两个速度根 \(\beta^\pm(r)\)；参考实现固定 \(\beta\)，求共振发生的位置 \(r(\beta)\)。

**代码里的 \(Q\) 和 \(1/(LQ)\)，都来自对半径积分的雅可比。**

先求 \(x(r)\) 的对数导数：
\[
\ln x(r)
=
\text{常数}-q\ln r
+\frac12\ln\left(1-\frac{r_s}{r}\right),
\]
所以
\[
\frac{d\ln x}{dr}
=
-\frac qr+\frac{r_s}{2r(r-r_s)}
=
-\frac{Q(r)}r,
\]
其中
\[
\boxed{
Q(r)=q-\frac{r_s}{2(r-r_s)}.
}
\tag{3}
\]
也就是
\[
\frac{dx}{dr}=-\frac{xQ}{r}.
\]

代码要求 \(Q>0\)，保证 \(x(r)\) 严格递减，因此每个速度至多对应一个共振半径。

接下来处理式 (1) 的 \(\delta\) 函数。因为
\[
\omega_D=\omega\frac{x}{g},
\]
所以在固定 \(r\) 或固定 \(\beta\) 的积分中，都可利用共振点上的雅可比写成
\[
\omega_D\delta(\omega-\omega_D)
=
x\,\delta(g-x).
\tag{4}
\]

将式 (4) 代入式 (1)，交换积分顺序：
\[
\tau
=
C\int d\beta\,
\frac{f(\beta)}{|\bar\beta|}
(1-\beta\mu)\mathcal P(\beta)
\int_{r_{\rm lo}}^{r_{\rm hi}}
\frac{dr}{rL(r)}\,x(r)\,
\delta\!\bigl(g(\beta)-x(r)\bigr).
\]

若共振半径 \(r_\beta=r(\beta)\) 位于传播区间内，则
\[
\delta\!\bigl(g(\beta)-x(r)\bigr)
=
\frac{\delta(r-r_\beta)}{|x'(r_\beta)|}.
\]
因此里面的半径积分恰好等于
\[
\frac{x(r_\beta)}{r_\beta L(r_\beta)}
\frac{1}{|x'(r_\beta)|}
=
\frac{x}{rL}\frac{r}{xQ}
=
\boxed{\frac1{L(r_\beta)Q(r_\beta)}}.
\]

于是得到文件注释中的第一条公式：
\[
\boxed{
\tau
=
C\int_{\mathcal S}d\beta\,
\frac{f(\beta)}{|\bar\beta|}
\frac{\mathcal P(\beta)(1-\beta\mu)}
{L[r(\beta)]Q[r(\beta)]}.
}
\tag{5}
\]

其中 \(\mathcal S\) 只包含分布支持范围内、且共振位置确实在传播区间内的速度：
\[
x(r_{\rm hi})\le g(\beta)\le x(r_{\rm lo}).
\tag{6}
\]

直观地说，\(1/L\) 来自局域路程与坐标半径的换算，\(1/Q\) 来自共振条件随半径变化的快慢。

**为什么这样就避开了两个速度根合并时的奇异性？**

如果先积掉速度，\(\delta\) 函数会引入
\[
\frac1{|g'(\beta)|},
\qquad
g'(\beta)=\gamma^3(\beta-\mu).
\]
在 \(\beta=\mu\) 时，这个导数为零，两个速度根合并，普通的根求和表达式出现奇异分母。

参考实现先积掉半径，使用的是
\[
\frac1{|x'(r)|}=\frac r{xQ}.
\]
只要 \(Q>0\)，它就没有这种奇异性。两个不同速度即使在同一个半径共振，也会在速度积分中分别计入，不需要显式求解 \(\beta^\pm\)。

**接下来推导为什么换成 \(u\) 后，贝塞尔函数消失了。**

项目使用单向 Boltzmann 分布。令
\[
s=\operatorname{sgn}(\beta_0),\qquad
a=\frac1{\gamma_0-1},\qquad
\gamma_0=\frac1{\sqrt{1-\beta_0^2}},
\]
在 \(0<s\beta<1\) 内，
\[
f(\beta)=\frac{\gamma^3e^{-a\gamma}}{K_1(a)}.
\tag{7}
\]

先计算平均速度的绝对值。利用
\[
d\gamma=|\beta|\gamma^3\,d|\beta|,
\]
得到
\[
|\bar\beta|
=
\frac1{K_1(a)}
\int_1^\infty e^{-a\gamma}\,d\gamma
=
\frac{e^{-a}}{aK_1(a)}.
\]
因此
\[
\boxed{
\frac{f(\beta)}{|\bar\beta|}
=
a\gamma^3e^{-a(\gamma-1)}.
}
\tag{8}
\]
这里的 \(K_1(a)\) 已经完全消掉了，所以参考实现不需要调用贝塞尔函数归一化。

再引入
\[
\boxed{u=\sqrt{a(\gamma-1)}},
\qquad
\gamma=1+\frac{u^2}{a}.
\]
由 \(|\beta|=\sqrt{1-\gamma^{-2}}\)，可得
\[
\boxed{
\beta(u)
=
s\,\frac{u\sqrt{2a+u^2}}{a+u^2}.
}
\tag{9}
\]
这就是代码中的 `beta(u)`。

由
\[
2u\,du=a\,d\gamma
=a|\beta|\gamma^3|d\beta|,
\]
得到雅可比
\[
\left|\frac{d\beta}{du}\right|
=
\frac{2u}{a|\beta|\gamma^3}.
\]
与式 (8) 相乘：
\[
\begin{aligned}
\frac{f(\beta)}{|\bar\beta|}
\left|\frac{d\beta}{du}\right|
&=
a\gamma^3e^{-u^2}
\frac{2u}{a|\beta|\gamma^3}\\
&=
\frac{2u}{|\beta|}e^{-u^2}\\
&=
\boxed{
\frac{2(a+u^2)}{\sqrt{2a+u^2}}e^{-u^2}.
}
\end{aligned}
\tag{10}
\]
这正是代码中的 `density`。它包含了 \(1/|\bar\beta|\)，**并不是单独归一化到 1 的概率密度**。

低温对应 \(a\gg1\)。此时
\[
u\approx\sqrt{\frac a2}\,|\beta|.
\]
速度分布虽然集中在很窄的 \(|\beta|\sim a^{-1/2}\) 范围内，变换后仍集中在 \(u\sim1\) 附近，因此容易被数值积分采样到。

**最终代码实际计算的就是下面这个积分。**

粒子静止系中的方向余弦为
\[
\mu'=\frac{\mu-\beta}{1-\beta\mu},
\]
偏振因子为
\[
\mathcal P(\beta)=
\begin{cases}
\dfrac12,&\text{E 模},\\[4pt]
\dfrac{\mu'^2}{2},&\text{O 模}.
\end{cases}
\]
于是
\[
\boxed{
\tau
=
C\int_{\mathcal U}du\,
\underbrace{\frac{2(a+u^2)}{\sqrt{2a+u^2}}e^{-u^2}}_{\texttt{density}}
\frac{
\underbrace{\mathcal P[\beta(u)]}_{\texttt{overlap}}\,
\underbrace{[1-\beta(u)\mu]}_{\texttt{doppler}}
}{
\underbrace{L[r(u)]}_{\texttt{lapse}}\,
Q[r(u)]
}.
}
\tag{11}
\]
这里 \(r(u)\) 通过
\[
x[r(u)]=g[\beta(u)]
\]
确定，\(\mathcal U\) 则由式 (6) 的传播区间限制确定。对应实现就是 [第 163–179 行](./test/resonance_reference.hpp:163)。

剩下的代码是在可靠地计算式 (11)：

- `log_x`、`log_g` 用对数计算共振条件，并通过二分法求 \(r(u)\)。
- `turning` 在 \(g\) 的极值点 \(\beta=\mu\) 分段。对应
  \[
  u_{\rm turn}
  =
  \sqrt{a\left(\frac1{\sqrt{1-\mu^2}}-1\right)}
  =
  \sqrt{\frac{a\mu^2}{\sqrt{1-\mu^2}\,[1+\sqrt{1-\mu^2}]}}.
  \]
- `breaks` 再找出共振半径等于两个传播端点的位置，准确划出 \(\mathcal U\)，避免漏掉很窄的有效区间。
- 在有效区间上做自适应 Gauss–Kronrod 积分；实际截断到 \(u=8\)，并把被截去的尾部误差上界计入 `result.error`。

所以它作为参考实现的独立性在于：**直接从含 \(\delta\) 函数的积分出发，解析积掉半径，再对变换后的粒子分布求积**；不依赖生产实现的速度根求解和 \(d\tau/dl\) 内核。
