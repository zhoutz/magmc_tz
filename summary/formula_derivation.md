# FD11 偏振演化的完整公式推导

**当前启动判据：** 已依照 `doc/adiabatic.md` 改为 $|d\chi_B/dl|>10^{-3}\kappa$，
并接入测地线 stepper 的事件。下文第 6.2 节保留 FD11 式 (34) 的历史推导；
当前实现与新增推导见 [fd11_adiabatic_event.md](fd11_adiabatic_event.md)。Stokes 表示和冻结判据继续沿用。

本文详细记录上一任务中 `src/dev/fd11.cpp` 所使用的公式推导，并逐项说明它们怎样对应到代码。主要依据为仓库中的 [FD11.pdf](../article/FD11.pdf)：Fernández & Davis (2011), *The X-Ray Polarization Signature of Quiescent Magnetars: Effect of Magnetospheric Scattering and Vacuum Polarization*, ApJ 730:131。

原任务的实现和验证记录见 [fd11.md](fd11.md)。本文侧重推导，不把已验证的数值结果当作数学证明，也不把程序采用的数值保护条件误称为论文原公式。

推导的逻辑顺序是：局部几何与单位 → 真空介电张量 → 电场传播矩阵 → Stokes 方程 → 启动与初值 → 基底变换 → 散射概率与光学深度 → 旋转积分器 → 冻结判据 → 径向解析解 → 输出与验证。

**阅读目录**

- [1. 记号、单位及近似范围](#section-1)
- [2. 路程、引力红移与偏振屏幕](#section-2)
- [3. 从 FD11 介电张量到横向传播矩阵](#section-3)
- [4. 本征模、双折射率及代码中的系数](#section-4)
- [5. 从两个复振幅逐项推导三个 Stokes 方程](#section-5)
- [6. 为什么延迟开始，以及如何设置初值](#section-6)
- [7. 屏幕坐标旋转与随磁场转动的方程](#section-7)
- [8. 从 FD11 式 (33) 推导偏振散射重叠](#section-8)
- [9. 共振速度根、原散射核及偏振相关光学深度](#section-9)
- [10. 用分部积分保留快速振荡的散射贡献](#section-10)
- [11. 散射事件位置、出射模式及频移](#section-11)
- [12. 从 Stokes 旋转到实四元数传播算符](#section-12)
- [13. 四阶 Magnus 步、系数来源和误差控制](#section-13)
- [14. 仅保存 Stokes 时，如何实现论文的复振幅冻结判据](#section-14)
- [15. 径向纯模光子的解析解与快速路径](#section-15)
- [16. 系数插值、求导与稠密偏振查询](#section-16)
- [17. 逃逸时的统一屏幕与输出](#section-17)
- [18. 与上述推导有关的原有几何公式补充](#section-18)
- [19. 独立复振幅验证方程与实现对应表](#section-19)

<a id="section-1"></a>

## 1. 记号、单位及近似范围

### 1.1 统一记号

| 记号 | 含义 |
|---|---|
| $s$ | 光子沿轨道走过的当地物理路程，单位 km |
| $r$ | Schwarzschild 面积半径，单位 km |
| $L(r)=\sqrt{1-r_s/r}$ | 引力红移因子 |
| $E_\infty$ | 无穷远光子能量，代码变量名 `omega_inf`，实际单位 keV |
| $E_{\rm loc}$ | 当地静止观测者测量的光子能量 |
| $k_0=E_{\rm loc}/(\hbar c)$ | 当地真空波数，单位 km$^{-1}$ |
| $A_x,A_y$ | 用于推导的横向复电场振幅，正式程序不把它们作为光子状态 |
| $S=(Q,U,V)$ | 每个光子的三个归一化实 Stokes 参数 |
| $\boldsymbol b=\boldsymbol B/B$ | 磁场单位方向 |
| $\phi_B$ | 横向磁场相对偏振屏幕 x 轴的角度 |
| $\Theta=2\phi_B$ | Stokes 平面中对应的双倍角，避免与磁余纬混淆 |
| $\kappa=k_0\Delta n$ | 两个本征模的相对相位增长率，单位 km$^{-1}$ |
| $h_0,d,g$ | 传播矩阵的迹部分、对角差的一半、非对角元 |
| $\xi=\omega_c/\omega_{\rm loc}$ | 共振频率比，对应代码 `Geometry::x` |
| $\mu=\hat k\cdot\hat B$ | 恒星静止系中光子与磁场夹角的余弦 |
| $\nu=(\mu-\beta)/(1-\beta\mu)$ | 散射粒子静止系中的入射方向余弦 |

全文用 $d,g$ 表示传播矩阵的两个系数，不复用论文介电系数 $a$ 或磁场单位向量 $\boldsymbol b$。速度 $\beta$ 以光速为单位。

### 1.2 本文推导采用的物理近似

1. 与 FD11 §2.1.3 一样，介电性质由真空极化主导，忽略等离子体的介电贡献。
2. 在实际积分偏振的区域使用弱场展开 $B\ll B_{\rm QED}$，只保留真空修正的一阶项。
3. 使用几何光学近似：背景变化尺度远大于波长，复包络变化慢于真空载波。
4. 在此精度下保留横向的 $A_x,A_y$，忽略纵向振幅对横向方程的二阶反馈。
5. 散射的速度分布、无反冲的频移及出射模抽样沿用原程序。新增工作是偏振的连续演化，以及耦合状态下的入射偏振重叠。

因此，“从复振幅到 Stokes 的等价”是指在上述传播矩阵下的严格代数等价；不是声称实现了任意强磁场、任意等离子体介电张量的完整方程。

<a id="section-2"></a>

## 2. 路程、引力红移与偏振屏幕

### 2.1 为什么 `omega_inf` 要除以红移因子

Schwarzschild 度规写为

\[
ds_{\rm spacetime}^2=-L^2c^2dt^2+L^{-2}dr^2+r^2d\Omega^2.
\]

时空静态性给出沿光路守恒的 $E_\infty$。当地静止观测者的固有时间满足 $d\tau=Ldt$，所以该观测者测量的频率比坐标时间对应的频率大 $1/L$：

\[
E_{\rm loc}=\frac{E_\infty}{L(r)}.
\]

波数由能量除以 $\hbar c$ 得到：

\[
\boxed{k_0=\frac{E_\infty}{L(r)\hbar c}.}
\]

代码的长度单位是 km，因而必须使用

\[
\hbar c=1.973269804\times10^{-13}\ {\rm keV\,km}.
\]

这里不能把以 keV 表示的 `omega_inf` 再误当作以 s$^{-1}$ 表示的角频率，否则会重复或遗漏一个 $\hbar$。

### 2.2 原程序的路程变量

轨道平面中的当地空间线元为

\[
d\ell^2=L^{-2}dr^2+r^2d\psi^2.
\]

程序以当地物理路程 $s=\ell$ 为独立变量。若传播方向与径向的夹角为 $\alpha$，单位方向的径向和切向分量分别是 $\cos\alpha,\sin\alpha$，因此

\[
\frac{1}{L}\frac{dr}{ds}=\cos\alpha,
\qquad r\frac{d\psi}{ds}=\sin\alpha.
\]

即

\[
\boxed{\frac{dr}{ds}=L\cos\alpha,\qquad
\frac{d\psi}{ds}=\frac{\sin\alpha}{r}.}
\]

这说明偏振方程中的 $d/ds$ 与原轨道积分使用的是同一个当地路程参数；不需要额外把 $dr$ 当作传播长度。

### 2.3 屏幕基底的构造

设轨道平面法向为固定的单位向量 $\boldsymbol n$，轨道平面内参考基为 $\boldsymbol e_1,\boldsymbol e_2$。定义

\[
\hat r=\cos\psi\,\boldsymbol e_1+\sin\psi\,\boldsymbol e_2,
\qquad
\hat\psi=\boldsymbol n\times\hat r,
\]

\[
\hat k=\cos\alpha\,\hat r+\sin\alpha\,\hat\psi.
\]

取横向屏幕基

\[
\boxed{\boldsymbol e_x=\boldsymbol n\times\hat k,
\qquad \boldsymbol e_y=\boldsymbol n.}
\]

因为 $\boldsymbol n\cdot\hat k=0$，所以 $e_x,e_y$ 都与 $\hat k$ 垂直，而且长度均为 1。再由向量三重积

\[
(\boldsymbol n\times\hat k)\times\boldsymbol n
=\hat k(\boldsymbol n\cdot\boldsymbol n)
-\boldsymbol n(\hat k\cdot\boldsymbol n)=\hat k,
\]

可知 $(e_x,e_y,\hat k)$ 是右手正交基。

屏幕系数的意义是**当地正交标架**中的分量，并不是把弯曲时空中的全局坐标向量当作欧氏向量做普通平行输运。

### 2.4 为什么这里不额外加入一个引力导致的屏幕内旋转项

利用球对称性可把轨道转到赤道面。此时轨道法向对应四维标架向量 $e_{\hat\vartheta}=r^{-1}\partial_\vartheta$。沿赤道光路，

\[
\partial_r(r^{-1})+\Gamma^\vartheta_{r\vartheta}r^{-1}
=-r^{-2}+r^{-1}r^{-1}=0,
\]

同时赤道上 $\cot\vartheta=0$，所以沿方位运动也不产生把该法向混入屏幕内另一个方向的分量。于是轨道法向沿光路平行输运。

再对 $e_x\cdot n=0$ 求协变导数，有

\[
n\cdot\nabla_k e_x=-e_x\cdot\nabla_k n=0.
\]

因此屏幕内没有 $e_x\leftrightarrow e_y$ 的额外转动。沿传播方向或时间方向的规范分量不改变观测到的横向偏振。这就是代码用这组屏幕基描述无双折射时的偏振平行输运的理由。

<a id="section-3"></a>

## 3. 从 FD11 介电张量到横向传播矩阵

### 3.1 弱场真空张量

FD11 式 (18)、(19) 为

\[
\varepsilon_{ij}=(1+a_{\rm vac})\delta_{ij}+q b_i b_j,
\qquad
\bar\mu_{ij}=(1+a_{\rm vac})\delta_{ij}+m b_i b_j.
\]

这里 $\bar\mu$ 是**逆磁导率张量**。弱场展开给出

\[
a_{\rm vac}=-2\delta,\quad q=7\delta,\quad m=-4\delta,
\qquad
\delta=\frac{\alpha_{\rm em}}{45\pi}
\left(\frac{B}{B_{\rm QED}}\right)^2.
\]

所有这些系数均为实数。各向同性部分的 $a_{\rm vac}$ 与之后使用的 Stokes 参数没有关系。

### 3.2 纵向振幅为什么不进入本阶横向方程

论文式 (23) 给出

\[
A_z=-\frac{\varepsilon_{zx}A_x+\varepsilon_{zy}A_y}{\varepsilon_{zz}}.
\]

其中

\[
\varepsilon_{zx}=q b_zb_x=O(\delta),\quad
\varepsilon_{zy}=q b_zb_y=O(\delta),\quad
\varepsilon_{zz}=1+O(\delta).
\]

所以 $A_z=O(\delta)A_\perp$。若它通过一个非对角介电项再反馈到横向方程，会多乘一个 $O(\delta)$，成为 $O(\delta^2)A_\perp$。本实现只保留一阶，因而舍去这种反馈。

### 3.3 逐项展开式 (21)

把论文式 (21) 按 $A_x,A_y$ 的系数写开：

\[
\begin{aligned}
\frac{dA_x}{ds}
={ik_0\over2}
\left(\bar\mu_{yy}-{\bar\mu_{xy}\bar\mu_{yx}\over\bar\mu_{xx}}\right)^{-1}
\Bigg\{&\left[\varepsilon_{xx}-\bar\mu_{yy}
+{\bar\mu_{yx}\over\bar\mu_{xx}}(\varepsilon_{yx}+\bar\mu_{xy})\right]A_x\\
+&\left[\varepsilon_{xy}+\bar\mu_{yx}
+{\bar\mu_{yx}\over\bar\mu_{xx}}(\varepsilon_{yy}-\bar\mu_{xx})\right]A_y\Bigg\}.
\end{aligned}
\]

下面分别数每一项的阶数。

首先，两个非对角逆磁导率的乘积为

\[
\bar\mu_{xy}\bar\mu_{yx}=m^2b_x^2b_y^2=O(\delta^2).
\]

所以大括号前的逆因子为 $1+O(\delta)$。而大括号内每个保留下来的系数已经是 $O(\delta)$，两者的修正乘积是二阶。因此本阶可以把这个逆因子置为 1。

其次，对角差为

\[
\begin{aligned}
\varepsilon_{xx}-\bar\mu_{yy}
&=(1+a_{\rm vac}+q b_x^2)
 -(1+a_{\rm vac}+m b_y^2)\\
&=q b_x^2-m b_y^2.
\end{aligned}
\]

可以看到 $1+a_{\rm vac}$ 在这里明确相消。

非对角和为

\[
\varepsilon_{xy}+\bar\mu_{yx}
=q b_xb_y+m b_yb_x=(q+m)b_xb_y.
\]

剩余两项各自包含两个一阶因子的乘积，例如

\[
{\bar\mu_{yx}\over\bar\mu_{xx}}
(\varepsilon_{yx}+\bar\mu_{xy})
=O(\delta)O(\delta)=O(\delta^2),
\]

因此都舍去。最后得到

\[
\frac{dA_x}{ds}={ik_0\over2}
\left[(q b_x^2-m b_y^2)A_x+(q+m)b_xb_y A_y\right].
\]

### 3.4 式 (22) 的展开

交换上述推导中的 x、y，可得

\[
\frac{dA_y}{ds}={ik_0\over2}
\left[(q+m)b_xb_y A_x+(q b_y^2-m b_x^2)A_y\right].
\]

因此横向 Jones 向量 $A=(A_x,A_y)^T$ 满足

\[
\boxed{A'=iHA,\qquad
H={k_0\over2}
\begin{pmatrix}
q b_x^2-m b_y^2&(q+m)b_xb_y\\
(q+m)b_xb_y&q b_y^2-m b_x^2
\end{pmatrix}.}
\]

矩阵 $H$ 为实对称矩阵，因而也是 Hermitian 矩阵。这里的正号 $+iH$ 来自论文载波与时间因子的约定，后续的 $V$ 符号必须与它保持一致。

<a id="section-4"></a>

## 4. 本征模、双折射率及代码中的系数

### 4.1 横向磁场角度

令光子与磁场的夹角为 $\vartheta_{kB}$，并记

\[
b_\perp=\sin\vartheta_{kB},\quad
b_x=b_\perp\cos\phi_B,\quad b_y=b_\perp\sin\phi_B.
\]

设

\[
e_O=\begin{pmatrix}\cos\phi_B\\\sin\phi_B\end{pmatrix},
\qquad
e_E=\begin{pmatrix}-\sin\phi_B\\\cos\phi_B\end{pmatrix}.
\]

O 模的电场沿横向磁场方向；E 模与它垂直。直接相乘可把传播矩阵写成

\[
H={k_0 b_\perp^2\over2}
\left(q e_Oe_O^T-m e_Ee_E^T\right).
\]

例如这个表达式的 xx 元为

\[
{k_0b_\perp^2\over2}
(q\cos^2\phi_B-m\sin^2\phi_B),
\]

正好等于上一节的 $H_{xx}$；xy 元也给出 $(q+m)\cos\phi_B\sin\phi_B$。

由于 $e_O,e_E$ 正交归一，立即得到本征值

\[
\lambda_O={k_0q b_\perp^2\over2},\qquad
\lambda_E=-{k_0m b_\perp^2\over2}.
\]

这些本征值是相对于真空载波的相位增长率，不是完整的真空波数。

### 4.2 相对相位率 $\kappa$

两模的本征值差为

\[
\lambda_O-\lambda_E
={k_0(q+m)b_\perp^2\over2}=k_0\Delta n.
\]

定义

\[
\boxed{\kappa=k_0\Delta n,\qquad
\Delta n={(q+m)\sin^2\vartheta_{kB}\over2}.}
\]

代入 $q+m=3\delta$，有

\[
\kappa={3\over2}k_0\delta\sin^2\vartheta_{kB}.
\]

由于 $B^2\sin^2\vartheta_{kB}=B_x^2+B_y^2$，再代入 $\delta$ 和 $k_0$：

\[
\begin{aligned}
\kappa
&={3\over2}{E_{\rm loc}\over\hbar c}
 {\alpha_{\rm em}\over45\pi}
 {B^2\over B_{\rm QED}^2}\sin^2\vartheta_{kB}\\
&=\boxed{{\alpha_{\rm em}\over30\pi\hbar c B_{\rm QED}^2}
 {E_\infty\over L(r)}(B_x^2+B_y^2)}.
\end{aligned}
\]

量纲检查：$E_\infty/(\hbar c)$ 为 km$^{-1}$，磁场平方比和红移因子均无量纲，所以 $\kappa$ 为 km$^{-1}$。这正是 `pol_coeff()` 中的表达式。

### 4.3 分离公共相位与偏振相关部分

写成

\[
H=h_0\mathbf1+
\begin{pmatrix}d&g\\g&-d\end{pmatrix},
\]

则根据矩阵元素定义

\[
h_0={H_{xx}+H_{yy}\over2},\quad
d={H_{xx}-H_{yy}\over2},\quad g=H_{xy}.
\]

依次代入：

\[
h_0={k_0\over4}(q-m)(b_x^2+b_y^2)
={11\over4}k_0\delta b_\perp^2={11\over6}\kappa,
\]

\[
d={k_0\over4}(q+m)(b_x^2-b_y^2)
={\kappa\over2}\cos2\phi_B,
\]

\[
g={k_0\over2}(q+m)b_xb_y
={\kappa\over2}\sin2\phi_B.
\]

这里使用了

\[
\cos2\phi_B={B_x^2-B_y^2\over B_x^2+B_y^2},\qquad
\sin2\phi_B={2B_xB_y\over B_x^2+B_y^2}.
\]

代码直接用这两个比值，不需要先算 $\phi_B$ 再计算三角函数。当 $B_\perp=0$ 时，$\kappa=0$，两个横向本征模退化；代码给出确定性的参考方向，而不除以零。

<a id="section-5"></a>

## 5. 从两个复振幅逐项推导三个 Stokes 方程

### 5.1 Stokes 定义及公共相位不变性

FD11 式 (36)–(39) 的约定为

\[
I=|A_x|^2+|A_y|^2,\quad Q=|A_x|^2-|A_y|^2,
\]

\[
U=2\operatorname{Re}(A_xA_y^*),\qquad
\boxed{V=2\operatorname{Im}(A_xA_y^*)}.
\]

如果同时把两个振幅乘上 $e^{i\chi}$，则

\[
(e^{i\chi}A_x)(e^{i\chi}A_y)^*=A_xA_y^*,
\]

而各分量模平方也不变，所以四个 Stokes 参数都不依赖公共相位。

把传播方程写成分量形式：

\[
A_x'=i(h_0+d)A_x+igA_y,
\qquad
A_y'=igA_x+i(h_0-d)A_y.
\]

以下不省略共轭项，以便明确检查正负号。

### 5.2 $I'$ 与 $Q'$

首先

\[
\begin{aligned}
(|A_x|^2)'
&=A_x'A_x^*+A_x(A_x')^*\\
&=[i(h_0+d)A_x+igA_y]A_x^*
 +A_x[-i(h_0+d)A_x^*-igA_y^*]\\
&=ig(A_yA_x^*-A_xA_y^*).
\end{aligned}
\]

令 $C=A_xA_y^*=(U+iV)/2$。于是

\[
A_yA_x^*-A_xA_y^*=C^*-C=-iV,
\]

所以

\[
(|A_x|^2)'=ig(-iV)=gV.
\]

同理，

\[
\begin{aligned}
(|A_y|^2)'
&=[igA_x+i(h_0-d)A_y]A_y^*
 +A_y[-igA_x^*-i(h_0-d)A_y^*]\\
&=ig(C-C^*)=ig(iV)=-gV.
\end{aligned}
\]

两式相加、相减，分别得到

\[
\boxed{I'=0,\qquad Q'=2gV.}
\]

### 5.3 $U'$ 与 $V'$

对 $C=A_xA_y^*$ 求导：

\[
\begin{aligned}
C'
&=A_x'A_y^*+A_x(A_y')^*\\
&=[i(h_0+d)A_x+igA_y]A_y^*
 +A_x[-igA_x^*-i(h_0-d)A_y^*]\\
&=i[(h_0+d)-(h_0-d)]C
 +ig|A_y|^2-ig|A_x|^2\\
&=2idC-igQ.
\end{aligned}
\]

代入 $C=(U+iV)/2$：

\[
C'=id(U+iV)-igQ=-dV+i(dU-gQ).
\]

另一方面 $C'=(U'+iV')/2$。比较实部和虚部：

\[
\boxed{U'=-2dV,\qquad V'=2dU-2gQ.}
\]

合并得到正式实现的三个实数方程：

\[
\boxed{
\begin{aligned}
Q'&=\kappa\sin\Theta\,V,\\
U'&=-\kappa\cos\Theta\,V,\\
V'&=\kappa(\cos\Theta\,U-\sin\Theta\,Q).
\end{aligned}}
\]

公共相位系数 $h_0$ 在每个等式中都相消，而不是被人为设为零。

### 5.4 叉乘形式与守恒量

定义

\[
\boldsymbol\Omega=(\kappa\cos\Theta,\kappa\sin\Theta,0).
\]

展开叉乘：

\[
\boldsymbol\Omega\times S
=\begin{pmatrix}
\kappa\sin\Theta V\\
-\kappa\cos\Theta V\\
\kappa\cos\Theta U-\kappa\sin\Theta Q
\end{pmatrix}.
\]

因此

\[
\boxed{S'=\boldsymbol\Omega\times S.}
\]

再直接求偏振长度的导数：

\[
\begin{aligned}
{d\over ds}(Q^2+U^2+V^2)
&=2Q(2gV)+2U(-2dV)+2V(2dU-2gQ)\\
&=4gQV-4dUV+4dUV-4gQV=0.
\end{aligned}
\]

对单个归一化纯态，恒等式也可直接从振幅证明：

\[
\begin{aligned}
Q^2+U^2+V^2
&=(|A_x|^2-|A_y|^2)^2+4|A_xA_y^*|^2\\
&=(|A_x|^2+|A_y|^2)^2=I^2=1.
\end{aligned}
\]

这个等式适用于每个纯态光子。把不同偏振的光子做统计平均后，平均 Stokes 向量的长度可以小于平均 $I$，不能再要求 ensemble 的偏振度为 1。

### 5.5 一个不能省略的符号检查

让横向磁场沿屏幕 x 轴，即 $\Theta=0$，并令 $\kappa$ 为常数。此时

\[
Q'=0,\quad U'=-\kappa V,\quad V'=\kappa U.
\]

从 $Q(0)=0,U(0)=1,V(0)=0$ 出发，解为

\[
U(s)=\cos(\kappa s),\qquad V(s)=\sin(\kappa s).
\]

独立地，Jones 解满足 $A_x\propto e^{i\lambda_Os}$、$A_y\propto e^{i\lambda_Es}$，所以 $A_xA_y^*\propto e^{i\kappa s}$，同样给出正的 $V=\sin(\kappa s)$。这同时检查了 $+iH$、Stokes 的共轭顺序和圆偏振符号。

<a id="section-6"></a>

## 6. 为什么延迟开始，以及如何设置初值

### 6.1 从相对相位得到 $\ell_A$

本征模之间的相对相位增量为

\[
d(\varphi_O-\varphi_E)=(\lambda_O-\lambda_E)ds=\kappa ds.
\]

在系数近似常量的局部区域，让相对相位改变一个弧度所需的长度为

\[
\boxed{\ell_A={1\over\kappa}={1\over k_0\Delta n}.}
\]

这就是论文式 (25)。它是“一弧度”尺度，不是完整一周的 $2\pi/\kappa$。

### 6.2 把论文式 (34) 改写为求根形式（历史启动判据）

论文选择

\[
\eta_{\rm couple}={\ell_A\over r}=10^{-3}
\]

作为开启数值积分的过渡面。向外传播时一般 $B$ 减弱，$\kappa$ 下降，$\ell_A/r$ 增大。因此积分开始的条件写为

\[
{1\over\kappa r}\ge\eta_{\rm couple}.
\]

当 $\kappa>0$ 时，乘上正数 $\kappa r$，得到

\[
1\ge\eta_{\rm couple}\kappa r
\quad\Longleftrightarrow\quad
F(s)=\eta_{\rm couple}\kappa(s)r(s)-1\le0.
\]

原实现使用后一个形式，避免显式计算可能发散的 $1/\kappa$。若发射或散射后的初始位置已经满足 $F\le0$，就在该位置开始；否则扫描轨道段并在首次检测到的越界区间内求根。

该历史实现默认扫描每个轨道步四个子区间。当前已移除此启动扫描，改为由测地线事件检测
$G=|d\chi_B/dl|-\epsilon_{\rm ad}\kappa$ 的过零，默认 $\epsilon_{\rm ad}=10^{-3}$，
轨道步上限仍为 $0.1r$。这是有限分辨率事件检测，不是关于所有可能非单调窄穿越的数学保证。

### 6.3 O 模和 E 模对应的 Stokes 初值

O 模可写为

\[
A=e^{i\chi_0}\begin{pmatrix}\cos\phi_B\\\sin\phi_B\end{pmatrix}.
\]

代入定义：

\[
Q_O=\cos^2\phi_B-\sin^2\phi_B=\cos2\phi_B,
\]

\[
U_O=2\cos\phi_B\sin\phi_B=\sin2\phi_B,
\qquad V_O=0.
\]

E 模可写为

\[
A=e^{i\chi_0}\begin{pmatrix}-\sin\phi_B\\\cos\phi_B\end{pmatrix},
\]

所以

\[
Q_E=-\cos2\phi_B,\qquad U_E=-\sin2\phi_B,\qquad V_E=0.
\]

最后

\[
\boxed{S_O=(\cos\Theta,\sin\Theta,0),\qquad S_E=-S_O.}
\]

任意公共初相位 $\chi_0$ 都消失。因此不必像复振幅实现那样再抽取一个随机公共相位。

在 `Mode` 阶段，`stokes={0,0,0}` 只是“尚未实例化这三个数”的程序占位值，不能把它解读为一个物理上完全不偏振的光子；物理状态由 O/E 标签表示。

<a id="section-7"></a>

## 7. 屏幕坐标旋转与随磁场转动的方程

### 7.1 从 Jones 分量推导 Stokes 的双倍角变换

设新的屏幕基相对旧基转过角度 $\chi$：

\[
\tilde e_x=\cos\chi\,e_x+\sin\chi\,e_y,
\qquad
\tilde e_y=-\sin\chi\,e_x+\cos\chi\,e_y.
\]

同一个电场向量在新基中的分量为

\[
\tilde A_x=cA_x+s_\chi A_y,\qquad
\tilde A_y=-s_\chi A_x+cA_y,
\]

其中 $c=\cos\chi,s_\chi=\sin\chi$。逐项展开模平方：

\[
|\tilde A_x|^2=c^2|A_x|^2+s_\chi^2|A_y|^2
+cs_\chi(A_xA_y^*+A_x^*A_y),
\]

\[
|\tilde A_y|^2=s_\chi^2|A_x|^2+c^2|A_y|^2
-cs_\chi(A_xA_y^*+A_x^*A_y).
\]

相减并用 $U=A_xA_y^*+A_x^*A_y$，得到

\[
\tilde Q=(c^2-s_\chi^2)Q+2cs_\chi U
=Q\cos2\chi+U\sin2\chi.
\]

再展开交叉乘积：

\[
\begin{aligned}
\tilde A_x\tilde A_y^*
&=-cs_\chi|A_x|^2+c^2A_xA_y^*
-s_\chi^2A_yA_x^*+cs_\chi|A_y|^2\\
&=-cs_\chi Q+(c^2-s_\chi^2){U\over2}
+i(c^2+s_\chi^2){V\over2}.
\end{aligned}
\]

分别取两倍实部和两倍虚部：

\[
\boxed{
\tilde Q=Q\cos2\chi+U\sin2\chi,\quad
\tilde U=-Q\sin2\chi+U\cos2\chi,\quad
\tilde V=V.}
\]

所以 Stokes 中出现双倍角，圆偏振在保持右手性的屏幕旋转下不变。这是坐标的被动变换；它与把偏振向量主动旋转相同角度的矩阵符号相反。

### 7.2 随磁场转动的基底

取 $\chi=\phi_B(s)$，新 x 轴沿磁场在屏幕内的投影，则

\[
Q_B=Q\cos\Theta+U\sin\Theta,
\qquad
U_B=-Q\sin\Theta+U\cos\Theta.
\]

必须对系数本身也求导。首先

\[
\begin{aligned}
Q_B'
&=Q'\cos\Theta+U'\sin\Theta
 +\Theta'(-Q\sin\Theta+U\cos\Theta)\\
&=\kappa\sin\Theta V\cos\Theta
 -\kappa\cos\Theta V\sin\Theta+\Theta'U_B\\
&=\Theta'U_B.
\end{aligned}
\]

其次

\[
\begin{aligned}
U_B'
&=-Q'\sin\Theta+U'\cos\Theta
 -\Theta'(Q\cos\Theta+U\sin\Theta)\\
&=-\kappa(\sin^2\Theta+\cos^2\Theta)V-\Theta'Q_B\\
&=-\kappa V-\Theta'Q_B.
\end{aligned}
\]

最后，由原方程

\[
V'=\kappa(-Q\sin\Theta+U\cos\Theta)=\kappa U_B.
\]

所以

\[
\boxed{
Q_B'=\Theta'U_B,\quad
U_B'=-\Theta'Q_B-\kappa V,\quad
V'=\kappa U_B.}
\]

写成叉乘形式：

\[
\boxed{S_B'=\boldsymbol\Omega_B\times S_B,
\qquad\boldsymbol\Omega_B=(\kappa,0,-\Theta').}
\]

这组方程仍然与原 Stokes 方程等价，没有删除模式耦合。耦合项就是 $\Theta'$。当 $\kappa\gg|\Theta'|$ 时，旋转轴接近固定的 Q 轴，因此更适合指数型数值积分。

### 7.3 不显式求角度也能计算 $\Theta'$

定义 $v_x=\kappa\cos\Theta,v_y=\kappa\sin\Theta$。求导：

\[
v_x'=\kappa'\cos\Theta-\kappa\Theta'\sin\Theta,
\]

\[
v_y'=\kappa'\sin\Theta+\kappa\Theta'\cos\Theta.
\]

计算交叉组合：

\[
\begin{aligned}
v_xv_y'-v_yv_x'
&=\kappa\kappa'\cos\Theta\sin\Theta
 +\kappa^2\Theta'\cos^2\Theta\\
&\quad-\kappa\kappa'\sin\Theta\cos\Theta
 +\kappa^2\Theta'\sin^2\Theta\\
&=\kappa^2\Theta'.
\end{aligned}
\]

又有 $v_x^2+v_y^2=\kappa^2$，故

\[
\boxed{\Theta'={v_xv_y'-v_yv_x'\over v_x^2+v_y^2}.}
\]

当前代码已移除系数多项式插值，统一使用局部角差的有限差分计算此导数。两个位置的角差通过

\[
\Delta\Theta=\operatorname{atan2}
(c_as_b-s_ac_b,\;c_ac_b+s_as_b)
\]

求得，其中 $c_a=\cos\Theta(a),s_a=\sin\Theta(a)$，这些系数均由 `pol_coeff` 直接计算。
然后用 $\Theta'\approx\Delta\Theta/(b-a)$；差分半宽为 `1e-4*stepper.h_old`，
端点裁剪到当前轨道稠密输出的定义域。这样不会直接对有分支跳跃的 `atan2` 输出做减法。

<a id="section-8"></a>

## 8. 从 FD11 式 (33) 推导偏振散射重叠

### 8.1 为什么使用粒子静止系的方向余弦

粒子沿磁场以 $\beta c$ 运动。光子在恒星静止系中平行磁场的动量为 $p_\parallel=E_{\rm loc}\mu/c$。沿磁场做 Lorentz 变换：

\[
E'=\gamma E_{\rm loc}(1-\beta\mu),\qquad
cp_\parallel'=\gamma E_{\rm loc}(\mu-\beta).
\]

相除得到

\[
\boxed{\nu={cp_\parallel'\over E'}
={\mu-\beta\over1-\beta\mu}.}
\]

两个共振速度根对应不同的 $\beta$，所以此处的 $\nu$ 需要分别计算。

### 8.2 为什么 O/E 基中的归一化偏振可以直接用于这个变换

下面用电磁波 Lorentz 变换说明，沿 B 的 boost 不会把 O、E 两个线偏振方向相互混合。

暂取 B 沿 Z 轴，光子方向为 $\hat k=(\sin\vartheta,0,\cos\vartheta)$，并选

\[
e_O=(\cos\vartheta,0,-\sin\vartheta),\qquad e_E=(0,1,0).
\]

这是一组右手横向基；若要让 O 轴严格指向 B 的投影，可以同时把这两个基向量取负，Stokes 不受影响。

令电场为 $\mathcal E=u e_O+v e_E$，并暂用 $c=1$。波磁场为

\[
\mathcal B=\hat k\times\mathcal E=u e_E-v e_O.
\]

沿 Z 方向 boost 时，

\[
\mathcal E_\perp'=\gamma(\mathcal E+\boldsymbol\beta\times\mathcal B)_\perp,
\qquad \mathcal E_Z'=\mathcal E_Z.
\]

逐分量得到

\[
\mathcal E_X'=\gamma u(\mu-\beta),\quad
\mathcal E_Y'=\gamma v(1-\beta\mu),\quad
\mathcal E_Z'=-u\sin\vartheta.
\]

设 Doppler 因子 $\mathcal D=\gamma(1-\beta\mu)$。光行差给出

\[
\cos\vartheta'=\nu,\qquad \sin\vartheta'={\sin\vartheta\over\mathcal D}.
\]

因此上面的三分量正好可以写成

\[
\mathcal E'=\mathcal D\left[u e_O'+v e_E'\right].
\]

两个 Jones 分量都乘上同一个实因子 $\mathcal D$，归一化到 $I=1$ 后，该因子消失。因而可以保留 O/E 基中的归一化 Stokes，只把方向余弦改成 $\nu$。

### 8.3 展开复振幅模平方

以磁场屏幕基为参考，FD11 式 (33) 可写为

\[
W_\sigma=|e_\sigma|^2={1\over2}|\nu A_{Bx}+i\sigma A_{By}|^2,
\qquad \sigma=\pm1.
\]

当前代码按论文对应的电荷符号约定取电子 $\sigma=-1$、正电子 $\sigma=+1$。下面的代数只依赖上式对 $\sigma$ 的定义。

令

\[
A_{Bx}=u+iv,\qquad A_{By}=w+iz.
\]

则

\[
\nu A_{Bx}+i\sigma A_{By}
=(\nu u-\sigma z)+i(\nu v+\sigma w).
\]

所以

\[
\begin{aligned}
2W_\sigma
&=(\nu u-\sigma z)^2+(\nu v+\sigma w)^2\\
&=\nu^2(u^2+v^2)+(w^2+z^2)+2\sigma\nu(vw-uz).
\end{aligned}
\]

另一方面

\[
V=2\operatorname{Im}[(u+iv)(w-iz)]=2(vw-uz).
\]

因此

\[
W_\sigma={1\over2}\left[\nu^2|A_{Bx}|^2+|A_{By}|^2\right]
+{\sigma\nu V\over2}.
\]

利用 $I=1$ 与 $Q_B=|A_{Bx}|^2-|A_{By}|^2$，解得

\[
|A_{Bx}|^2={1+Q_B\over2},\qquad
|A_{By}|^2={1-Q_B\over2}.
\]

代回并整理：

\[
\boxed{W_\sigma={1\over4}
\left[(1+\nu^2)+(\nu^2-1)Q_B\right]
+{\sigma\nu V\over2}.}
\]

这个推导说明圆偏振项的系数是 $\sigma\nu/2$，不是 $\sigma\nu$，其符号也不能脱离前面的 $V=+2\operatorname{Im}(A_xA_y^*)$ 约定。

### 8.4 纯模极限与正性检查

O 模 $Q_B=1,V=0$，所以

\[
W_O={1\over4}[(1+\nu^2)+(\nu^2-1)]={\nu^2\over2}.
\]

E 模 $Q_B=-1,V=0$，所以

\[
W_E={1\over4}[(1+\nu^2)-(\nu^2-1)]={1\over2}.
\]

由其原始模平方形式显然 $W_\sigma\ge0$。再用 Cauchy–Schwarz 不等式，

\[
|\nu A_{Bx}+i\sigma A_{By}|^2
\le(\nu^2+1)(|A_{Bx}|^2+|A_{By}|^2)\le2,
\]

所以 $W_\sigma\le1$。代码的 `[0,1]` 截断是抑制浮点误差越界，不是额外的物理修正。

<a id="section-9"></a>

## 9. 共振速度根、原散射核及偏振相关光学深度

### 9.1 把共振条件化成二次方程

电子静止系中的共振条件为

\[
\gamma\omega_{\rm loc}(1-\beta\mu)=\omega_c.
\]

令 $\xi=\omega_c/\omega_{\rm loc}>0$，移项得

\[
1-\beta\mu=\xi\sqrt{1-\beta^2}.
\]

两边平方：

\[
1-2\beta\mu+\beta^2\mu^2=\xi^2-\xi^2\beta^2,
\]

整理为

\[
\boxed{(\xi^2+\mu^2)\beta^2-2\mu\beta+1-\xi^2=0.}
\]

判别式的一部分为

\[
\begin{aligned}
\mu^2-(\xi^2+\mu^2)(1-\xi^2)
&=\mu^2-\xi^2-\mu^2+\xi^4+\xi^2\mu^2\\
&=\xi^2(\xi^2+\mu^2-1).
\end{aligned}
\]

定义 $D=\xi^2+\mu^2-1$，则两个形式上的根为

\[
\boxed{\beta_\pm={\mu\pm\xi\sqrt D\over\xi^2+\mu^2}.}
\]

只有 $D\ge0$，并且根落在粒子分布支集内时，才有相应散射贡献。代码用稳定的二次方程求根器，而不是直接依赖可能发生抵消的分子表达式。

### 9.2 为什么原 O 模因子可以写成 $D/(2\xi^2)$

根据光行差公式，

\[
1-\nu^2
={(1-\beta\mu)^2-(\mu-\beta)^2\over(1-\beta\mu)^2}.
\]

展开分子：

\[
1-2\beta\mu+\beta^2\mu^2-\mu^2+2\beta\mu-\beta^2
=(1-\mu^2)(1-\beta^2).
\]

在共振上 $(1-\beta\mu)^2=\xi^2(1-\beta^2)$，所以

\[
1-\nu^2={1-\mu^2\over\xi^2},\qquad
\boxed{\nu^2={D\over\xi^2}.}
\]

这个平方对两个速度根相同。因此纯 O 模有 $W_O=D/(2\xi^2)$，纯 E 模有 $W_E=1/2$，可以提到两个根的和外面。但耦合后的 $\sigma\nu V/2$ 依赖 $\nu$ 的符号，不能再这样提取。

### 9.3 原程序中速度权重 $f(\beta)(1-\beta^2)^{3/2}$ 的来历

把粒子数密度写为 $n_e$，速度概率密度写为 $f(\beta)$。从论文共振截面中的结构

\[
\sigma_{\rm res}={4\pi^2e\over B}(1-\beta\mu)W_\sigma\,
\omega_D\,\delta(\omega_{\rm loc}-\omega_D),
\]

出发，其中

\[
\omega_D={\omega_c\sqrt{1-\beta^2}\over1-\beta\mu}.
\]

对其对数求导：

\[
\begin{aligned}
{1\over\omega_D}{d\omega_D\over d\beta}
&=-{\beta\over1-\beta^2}+{\mu\over1-\beta\mu}\\
&={-\beta(1-\beta\mu)+\mu(1-\beta^2)
\over(1-\beta^2)(1-\beta\mu)}\\
&={\mu-\beta\over(1-\beta^2)(1-\beta\mu)}.
\end{aligned}
\]

利用 delta 函数换元公式，速度积分变为各根的求和：

\[
\int F(\beta)\delta(\omega_{\rm loc}-\omega_D)d\beta
=\sum_i{F(\beta_i)\over|\omega_D'(\beta_i)|}.
\]

因此每个根的运动学因子是

\[
{(1-\beta\mu)\omega_D\over|\omega_D'|}
={(1-\beta^2)(1-\beta\mu)^2\over|\mu-\beta|}.
\]

在共振上，

\[
1-\beta\mu=\xi\sqrt{1-\beta^2},
\]

且由上一小节

\[
|\mu-\beta|=|\nu|(1-\beta\mu)
={\sqrt D\over\xi}\,\xi\sqrt{1-\beta^2}
=\sqrt D\sqrt{1-\beta^2}.
\]

所以

\[
{(1-\beta^2)(1-\beta\mu)^2\over|\mu-\beta|}
={\xi^2(1-\beta^2)^2\over\sqrt D\sqrt{1-\beta^2}}
={\xi^2\over\sqrt D}(1-\beta^2)^{3/2}.
\]

这给出粒子分布和共振 Jacobian 合并后的权重

\[
R_i=f(\beta_i)(1-\beta_i^2)^{3/2}.
\]

当前代码已删除独立的 `rate()` 函数。`polarized_rates()` 在同一循环中计算此权重，
乘上当前偏振状态的重叠因子后返回 $R_iW_{\sigma,i}$。

### 9.4 空间前因子与完整光学深度

原模型中单一电荷组分承担电流，粒子密度与扭转磁场满足

\[
{e n_e\over B}={p+1\over4\pi r|\bar\beta|}
{B_\varphi\over B_\vartheta}.
\]

因此 $n_e(4\pi^2 e/B)$ 等于

\[
{(p+1)\pi\over r|\bar\beta|}{B_\varphi\over B_\vartheta}.
\]

再乘上上一节的 $\xi^2/\sqrt D$，得到代码 `basepref`：

\[
\mathcal P(s)={(p+1)\pi\over r|\bar\beta|}
{B_\varphi\over B_\vartheta}{\xi^2\over\sqrt D}.
\]

故光学深度增长率为

\[
\boxed{{d\tau\over ds}=\mathcal P(s)\sum_i R_i(s)W_{\sigma,i}(s).}
\]

`Geometry` 只保存不含重叠因子的 `basepref`，不再另存 `pref`。
偏振 ODE 启动前或关闭偏振演化时，`polarized_rates()` 使用纯模式重叠
$W_E=1/2$、$W_O=D/(2\xi^2)$；启动后使用 Stokes 重叠表达式。
因此所有阶段的直接光学深度计算均为 `basepref * (rates[0] + rates[1])`。
`opacity_terms()` 为分部积分分别计算下面的三个线性系数，不提前乘入具体 Stokes 状态的重叠。

沿用 km 作为 r 的单位时，整个几何前因子就是 km$^{-1}$，与 ds 的单位匹配。

将重叠表达式代入，收集与 $1,Q_B,V$ 成比例的项：

\[
{d\tau\over ds}=c_0+c_QQ_B+c_VV,
\]

\[
\boxed{
\begin{aligned}
c_0&=\mathcal P\sum_iR_i{1+\nu_i^2\over4},\\
c_Q&=\mathcal P\sum_iR_i{\nu_i^2-1\over4},\\
c_V&=\mathcal P\sum_iR_i{\sigma\nu_i\over2}.
\end{aligned}}
\]

这些量就是 `opacity_terms()` 的三个输出；它们是实光学深度系数，不是电场振幅。

<a id="section-10"></a>

## 10. 用分部积分保留快速振荡的散射贡献

### 10.1 普通求积为何可能变慢

当某条光路在 $B\parallel k$ 附近开启模式耦合，随后进入较强横向磁场区域时，$\kappa$ 可能再次很大。此时 $U_B,V$ 快速振荡。即使单步的 Stokes 传播可以用矩阵指数高效完成，直接对 $c_VV$ 做普通求积仍可能需要解析许多周期。

本实现没有令 $V=0$，也没有把这个项直接做相位平均，而是利用它满足的微分方程改写积分。

### 10.2 分部积分的每一步

从已经推导的

\[
U_B'=-\Theta'Q_B-\kappa V
\]

解出

\[
V=-{U_B'\over\kappa}-{\Theta'Q_B\over\kappa}.
\]

乘以 $c_V$ 并在区间 $[l,r]$ 积分：

\[
\int_l^r c_VV\,ds
=-\int_l^r{c_V\over\kappa}U_B'\,ds
-\int_l^r{c_V\Theta'\over\kappa}Q_B\,ds.
\]

对第一个积分使用乘积求导恒等式

\[
\left({c_V\over\kappa}U_B\right)'
=\left({c_V\over\kappa}\right)'U_B
+{c_V\over\kappa}U_B',
\]

于是

\[
-\int_l^r{c_V\over\kappa}U_B'\,ds
=-\left[{c_V\over\kappa}U_B\right]_l^r
+\int_l^r\left({c_V\over\kappa}\right)'U_B\,ds.
\]

合并得到

\[
\boxed{
\int_l^r c_VV\,ds
=-\left[{c_V\over\kappa}U_B\right]_l^r
+\int_l^r\left[
\left({c_V\over\kappa}\right)'U_B
-{c_V\Theta'\over\kappa}Q_B\right]ds.}
\]

所以完整光学深度为

\[
\boxed{
\Delta\tau=\int_l^r\left[c_0+c_QQ_B
+\left({c_V\over\kappa}\right)'U_B
-{c_V\Theta'\over\kappa}Q_B\right]ds
+{c_V(l)U_B(l)\over\kappa(l)}
-{c_V(r)U_B(r)\over\kappa(r)}.}
\]

这是一个恒等变换。在系数变化慢于偏振相位的区域，新的振荡项因 $1/\kappa$ 而较小。实际程序仍有有限差分、插值和数值求积误差，不能把“恒等变换”理解为“整个数值计算没有误差”。

### 10.3 有限差分与不连续端点

代码用局部有限差分估计

\[
\left({c_V\over\kappa}\right)'(s)
\approx{{c_V(b)\over\kappa(b)}-{c_V(a)\over\kappa(a)}\over b-a},
\]

并用第 7.3 节的角差除以 $b-a$ 估计 $\Theta'$。采样点限制在同一个已切分积分区间的内部。

这一步的重要原因是原程序的单向 Boltzmann 分布：以 $\beta<0$ 的情形为例，支集为 $(-1,0)$，而 $\beta=0$ 及外部的密度定义为零。支集内部

\[
f(\beta)={e^{-a_T(\gamma-1)}\over K\,(1-\beta^2)^{3/2}},
\qquad \gamma={1\over\sqrt{1-\beta^2}},
\]

其中 $K$ 是代码中的正归一化常数。故

\[
\lim_{\beta\to0^-}f(\beta)={1\over K},\qquad f(0)=0.
\]

当某个共振根穿过 $\beta=0$，光学深度系数可能存在有限跳跃。如果直接取切点上由浮点舍入决定的那一侧系数，分部积分边界项就可能错一整个跳跃量。

正确做法是在每个平滑区间单独使用恒等式，左端取 $l^+$ 的系数、右端取 $r^-$ 的系数。代码用很小的内移量近似这些单侧极限，有限差分也不跨越支集切点。

### 10.4 实现中的启用条件

当前实现仅在已开始积分且尚未冻结的区间、估计相位 $\kappa(s_{\rm mid})(r-l)>256$、且两个端点的共振判别量 $D>10^{-8}$ 时使用分部积分。否则使用原始散射率求积。

256 rad 和 $10^{-8}$ 是数值算法的切换参数，不是 FD11 的新物理判据。支集边界、共振根合并点、偏振稠密解分段边界均被用于切分求积区间。

<a id="section-11"></a>

## 11. 散射事件位置、出射模式及频移

### 11.1 随机光学深度

若 $u$ 在 $(0,1)$ 上均匀分布，定义目标光学深度 $\tau_*=-\ln u$。因为

\[
P(\tau_*>t)=P(u<e^{-t})=e^{-t},
\]

它具有所需的指数生存概率。

沿光路累积 $\tau(s)=\int_{s_0}^s(d\tau/ds')ds'$，散射位置满足

\[
\boxed{\ln u+\int_{s_0}^{s_{\rm sc}}{d\tau\over ds'}ds'=0.}
\]

代码中初始 `tau=log(ran.U())`，再不断加上正的段光学深度，检测其何时过零。求根函数的导数就是原始散射率 $d\tau/ds$。

### 11.2 选择哪个共振速度根

在固定散射位置，两根共用同一空间前因子 $\mathcal P$，所以选择第 i 根的概率为

\[
P_i={R_iW_{\sigma,i}\over\sum_jR_jW_{\sigma,j}}.
\]

纯模区可以消去两根共有的重叠；耦合区必须保留各根不同的重叠，尤其是 $\sigma\nu_iV/2$ 项。

### 11.3 出射 O/E 模概率

沿用 FD11 与原代码的无反冲共振散射模型，固定粒子静止系的出射方向余弦 $\nu_{\rm out}$ 后，出射纯 O/E 模的重叠分别为

\[
W_{O,\rm out}={\nu_{\rm out}^2\over2},\qquad
W_{E,\rm out}={1\over2}.
\]

归一化后

\[
\boxed{P(E)={1\over1+\nu_{\rm out}^2},\qquad
P(O)={\nu_{\rm out}^2\over1+\nu_{\rm out}^2}.}
\]

所以每次散射后重新抽取一个纯模并重置为 `Mode` 状态；前一段连续演化的相干态不直接保留到新的出射方向上。这是所沿用的散射模型规定，不能从“无吸收的传播守恒”单独推出。

### 11.4 出射方向抽样为何对应 $1+\nu_{\rm out}^2$

对出射两模求和，方向概率密度正比于 $1+\nu_{\rm out}^2$。由于

\[
\int_{-1}^{1}(1+\nu^2)d\nu=2+{2\over3}={8\over3},
\]

归一化密度为

\[
p(\nu)={3\over8}(1+\nu^2).
\]

它可分解为：概率 $3/4$ 抽一个 $[-1,1]$ 均匀变量，概率 $1/4$ 抽密度 $3\nu^2/2$ 的变量。后一分布可由三个 $[0,1]$ 均匀变量的最大值再随机赋予正负号获得，因为

\[
P(\max(u_1,u_2,u_3)<m)=m^3,
\qquad p_{\max}(m)=3m^2.
\]

因此混合密度为

\[
{3\over4}{1\over2}+{1\over4}{3\nu^2\over2}
={3\over8}(1+\nu^2),
\]

对应原 `sample_mup()` 的抽样方式。

### 11.5 方向和能量变回恒星静止系

对光行差关系解出出射静止系余弦：

\[
\nu_{\rm out}={\mu_{\rm out}-\beta\over1-\beta\mu_{\rm out}}
\quad\Longrightarrow\quad
\boxed{\mu_{\rm out}={\nu_{\rm out}+\beta\over1+\beta\nu_{\rm out}}.}
\]

忽略反冲时，粒子静止系散射前后能量相同，因此

\[
\gamma E_{{\rm loc,in}}(1-\beta\mu_{\rm in})
=\gamma E_{{\rm loc,out}}(1-\beta\mu_{\rm out}).
\]

消去 $\gamma$，并注意散射发生在同一位置，入射和出射的红移因子相同，得到

\[
\boxed{E_{\infty,\rm out}=E_{\infty,\rm in}
{1-\beta\mu_{\rm in}\over1-\beta\mu_{\rm out}}.}
\]

这些关系保留了原代码的散射行为；偏振状态的新增内容主要进入入射散射率和之后的传播过程。

<a id="section-12"></a>

## 12. 从 Stokes 旋转到实四元数传播算符

### 12.1 常系数方程的精确解是一个旋转

若一小段内 $\boldsymbol\Omega$ 固定，令 $\omega=|\boldsymbol\Omega|$、$\hat u=\boldsymbol\Omega/\omega$。把 Stokes 分成

\[
S_\parallel=(\hat u\cdot S)\hat u,\qquad
S_\perp=S-S_\parallel.
\]

因为叉乘与 $\hat u$ 垂直，所以 $S_\parallel'=0$。对垂直部分再求一次导数：

\[
S_\perp''=\omega^2\hat u\times(\hat u\times S_\perp)
=\omega^2[\hat u(\hat u\cdot S_\perp)-S_\perp]
=-\omega^2S_\perp.
\]

所以经过长度 h 后

\[
S(s+h)=S_\parallel+\cos(\omega h)S_\perp
+\sin(\omega h)\hat u\times S_\perp.
\]

令 $\zeta=\omega h$，重新组合便得到 Rodrigues 公式：

\[
\boxed{S_{\rm new}=\cos\zeta\,S
+(1-\cos\zeta)\hat u(\hat u\cdot S)
+\sin\zeta\,\hat u\times S.}
\]

它对任意 $\zeta$ 都是一个长度守恒的旋转，不会因为相位大就出现显式 Runge–Kutta 方法那样的线性稳定性限制。但背景变化仍会带来截断误差，所以仍然需要自适应精度控制。

### 12.2 四元数形式以及代码中的叉乘更新

用单位四元数

\[
w=\cos{\zeta\over2},\qquad
\boldsymbol v=\hat u\sin{\zeta\over2}
\]

表示该旋转。显然 $w^2+|v|^2=1$。由倍角公式

\[
\cos\zeta=w^2-|v|^2,\quad
1-\cos\zeta=2|v|^2,\quad
\sin\zeta\,\hat u=2wv,
\]

Rodrigues 公式成为

\[
S_{\rm new}=(w^2-|v|^2)S+2v(v\cdot S)+2w(v\times S).
\]

再用 $w^2+|v|^2=1$ 与

\[
v\times(v\times S)=v(v\cdot S)-|v|^2S,
\]

得到

\[
\boxed{S_{\rm new}=S+2w(v\times S)+2v\times(v\times S).}
\]

代码先计算 $t=2v\times S$，再返回 $S+wt+v\times t$，正是上式。

对于输入的旋转向量 $a=h\boldsymbol\Omega$，需要计算

\[
v={\sin(|a|/2)\over|a|}a.
\]

当 $|a|$ 很小时，Taylor 展开为

\[
{\sin(|a|/2)\over|a|}
={1\over2}-{|a|^2\over48}+O(|a|^4),
\]

这解释了 `Rotation::exponential()` 中的短级数分支。

### 12.3 与 Jones 算符的对应关系

为稍后推导冻结指标，需要保留旋转的 Jones 对应，而不能只知道它在 Stokes 球上的几何效果。

定义通常的 Pauli 矩阵

\[
\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
\sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad
\sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
\]

由 Stokes 的定义，

\[
A^\dagger\sigma_zA=Q,\qquad
A^\dagger\sigma_xA=U,
\]

而

\[
A^\dagger\sigma_yA=-iA_x^*A_y+iA_y^*A_x=-V.
\]

因此适合本程序 $S=(Q,U,V)$ 排列的三个矩阵是

\[
T_Q=\sigma_z,\qquad T_U=\sigma_x,\qquad T_V=-\sigma_y,
\]

满足 $A^\dagger T_iA=S_i$。注意 $T_V$ 前面的负号不能省略。

它们的乘法规则为

\[
T_iT_j=\delta_{ij}\mathbf1-i\epsilon_{ijk}T_k.
\]

例如 $T_QT_U=\sigma_z\sigma_x=i\sigma_y=-iT_V$，说明此排列的符号与通常 $(\sigma_x,\sigma_y,\sigma_z)$ 的循环顺序不同。

相应的单位 Jones 算符是

\[
\boxed{G=w\mathbf1+i(v_QT_Q+v_UT_U+v_VT_V).}
\]

令 $K=v\cdot T$。上面的乘法规则给出 $K^2=|v|^2\mathbf1$，因为对称的 $v_iv_j$ 与反对称的 $\epsilon_{ijk}$ 收缩为零。因此

\[
G^\dagger G=(w\mathbf1-iK)(w\mathbf1+iK)
=(w^2+|v|^2)\mathbf1=\mathbf1.
\]

在一个无穷小步中，$v=h\boldsymbol\Omega/2+O(h^3)$。计算

\[
G^\dagger T_iG
=T_i+i[T_i,v_jT_j]+O(|v|^2)
=T_i+2\epsilon_{ijk}v_jT_k+O(|v|^2),
\]

故 $\Delta S=2v\times S=h\boldsymbol\Omega\times S+O(h^2)$。这检查了上述 Jones 算符对应的是正确方向的 Stokes 旋转。

### 12.4 连续两个传播步怎样合成

先应用 $G_1=w_1\mathbf1+iv_1\cdot T$，再应用 $G_2=w_2\mathbf1+iv_2\cdot T$，总算符为 $G_2G_1$。展开：

\[
G_2G_1=w_2w_1\mathbf1+i(w_2v_1+w_1v_2)\cdot T
-(v_2\cdot T)(v_1\cdot T).
\]

利用

\[
(v_2\cdot T)(v_1\cdot T)
=(v_2\cdot v_1)\mathbf1-i(v_2\times v_1)\cdot T,
\]

得到

\[
\boxed{
w_{21}=w_2w_1-v_2\cdot v_1,\qquad
v_{21}=w_2v_1+w_1v_2+v_2\times v_1.}
\]

这就是 `compose(after,before)` 的顺序。若另外保留公共相位 $e^{i\Phi_1},e^{i\Phi_2}$，由于它们是单位矩阵的标量倍数，直接相加 $\Phi_{21}=\Phi_2+\Phi_1$ 即可。

<a id="section-13"></a>

## 13. 四阶 Magnus 步、系数来源和误差控制

### 13.1 反对称矩阵生成元

把 $S'=\boldsymbol\Omega\times S$ 写成线性系统

\[
S'=\mathcal L(s)S,\qquad
\mathcal L(s)=
\begin{pmatrix}
0&-\Omega_V&\Omega_U\\
\Omega_V&0&-\Omega_Q\\
-\Omega_U&\Omega_Q&0
\end{pmatrix}.
\]

矩阵 $\mathcal L^T=-\mathcal L$，因此它的指数是正交旋转。背景变化时不同位置的 $\mathcal L$ 通常不对易，不能简单把所有系数平均一次就当作精确解。

### 13.2 中点展开与 Magnus 对易子项

考虑长度 h 的一步，以中点为零点，令局部坐标 $t\in[-h/2,h/2]$，并写

\[
\mathcal L(t)=D_0+tD_1+{t^2\over2}D_2+\cdots.
\]

传播算符的 Magnus 对数前两项为

\[
\mathcal M_1=\int_{-h/2}^{h/2}\mathcal L(t)dt,
\]

\[
\mathcal M_2={1\over2}\int_{-h/2}^{h/2}dt_1
\int_{-h/2}^{t_1}dt_2\,[\mathcal L(t_1),\mathcal L(t_2)].
\]

首先，奇函数项的对称积分为零，且

\[
\int_{-h/2}^{h/2}t^2dt={h^3\over12},
\]

因此

\[
\mathcal M_1=hD_0+{h^3\over24}D_2+O(h^5).
\]

对第二项，只保留最低的非零对易子：

\[
[D_0+t_1D_1,D_0+t_2D_1]=(t_2-t_1)[D_0,D_1].
\]

内层积分为

\[
\int_{-h/2}^{t_1}(t_2-t_1)dt_2=-{1\over2}(t_1+h/2)^2.
\]

再作外层积分并乘上前面的 $1/2$，得到

\[
\mathcal M_2=-{h^3\over12}[D_0,D_1]+O(h^5).
\]

还需要检查第三及更高的 Magnus 项是否会贡献同阶项。第三项包含三重积分和嵌套对易子；若三个生成元都取常数 $D_0$，嵌套对易子为零，因此非零项至少还要包含一个 $tD_1$，其阶数至少为 $h^4$。第四及更高项同理不会产生低于 $h^5$ 的新贡献。

在以同一个中点表达时，把 h 换成 $-h$ 相当于逆向传播，传播算符变成逆算符，因此其对数变号。这个对称性排除了四阶的偶数 h 项。于是四阶方法需要匹配

\[
\mathcal M=hD_0+{h^3\over24}D_2
-{h^3\over12}[D_0,D_1]+O(h^5).
\]

### 13.3 两个指数的系数怎样确定

取两个 Gauss 点 $t_\pm=\pm\gamma_G h$，其中 $\gamma_G=\sqrt3/6$。记 $\mathcal L_\pm=\mathcal L(t_\pm)$。考虑先后两个指数

\[
X=h(a_M\mathcal L_-+b_M\mathcal L_+),
\qquad
Y=h(b_M\mathcal L_-+a_M\mathcal L_+),
\]

总步算符为 $e^Ye^X$，右边的 X 先作用。

要匹配积分的一阶项，首先要求 $a_M+b_M=1/2$。于是

\[
X+Y={h\over2}(\mathcal L_-+\mathcal L_+).
\]

展开两采样点，

\[
\mathcal L_-+\mathcal L_+=2D_0+\gamma_G^2h^2D_2+O(h^4).
\]

因为 $\gamma_G^2=1/12$，得到

\[
X+Y=hD_0+{h^3\over24}D_2+O(h^5),
\]

与 $\mathcal M_1$ 一致。

接着使用 Baker–Campbell–Hausdorff 展开

\[
\log(e^Ye^X)=X+Y+{1\over2}[Y,X]+\cdots.
\]

直接计算

\[
[Y,X]=(b_M^2-a_M^2)h^2[\mathcal L_-,\mathcal L_+]
=-(a_M-b_M)(a_M+b_M)h^2[\mathcal L_-,\mathcal L_+].
\]

而

\[
[\mathcal L_-,\mathcal L_+]
=2\gamma_G h[D_0,D_1]+O(h^3).
\]

所以

\[
{1\over2}[Y,X]=-{(a_M-b_M)\gamma_G\over2}
h^3[D_0,D_1]+O(h^5).
\]

要使它等于 Magnus 的 $-h^3[D_0,D_1]/12$，必须有

\[
a_M-b_M={1\over6\gamma_G}={1\over\sqrt3}.
\]

联立和与差，得

\[
\boxed{a_M={3+2\sqrt3\over12},\qquad
b_M={3-2\sqrt3\over12}.}
\]

剩余嵌套对易子的最低组合为

\[
{1\over12}\left([Y,[Y,X]]+[X,[X,Y]]\right)
={1\over12}[Y-X,[Y,X]].
\]

因为 $Y-X=O(h^2)$、$[Y,X]=O(h^3)$，这个组合从 $O(h^5)$ 才开始。因此上述匹配给出局部误差 $O(h^5)$、整体四阶的方法。

### 13.4 实际应用到 Stokes 时的两个旋转向量

把矩阵线性组合重新写成叉乘生成元，只需构造

\[
u_1=h(a_M\boldsymbol\Omega_-+b_M\boldsymbol\Omega_+),
\qquad
u_2=h(b_M\boldsymbol\Omega_-+a_M\boldsymbol\Omega_+).
\]

先按 $u_1$ 旋转，再按 $u_2$ 旋转，并按第 12.4 节合成对应四元数。这样没有显式计算矩阵对易子，却保留了四阶所需的对易子效应。

公共相位独立用两点 Gauss 积分：

\[
\Phi=\int_s^{s+h}h_0(s')ds'
\approx{h\over2}(h_{0,-}+h_{0,+})
=\boxed{{11h\over12}(\kappa_-+\kappa_+)}.
\]

### 13.5 在转动基底中积分后，怎样换回原屏幕

令 $R_V(\Theta)$ 表示绕 Stokes 的 V 轴主动旋转角 $\Theta$ 的矩阵。第 7 节的被动变换就是

\[
S_B(s)=R_V(-\Theta(s))S(s).
\]

若转动基底中的数值传播算符为 $\mathcal U_B$，则

\[
\boxed{S(s+h)=R_V(\Theta_{\rm end})\mathcal U_B
R_V(-\Theta_{\rm begin})S(s).}
\]

程序把它改写为

\[
R_V(\Theta_{\rm begin})
\left[R_V(\Delta\Theta)\mathcal U_B\right]
R_V(-\Theta_{\rm begin}),
\]

先组合局部角差的旋转，再把四元数向量部转回起始屏幕。这样避免独立求两个半角时，`atan2` 分支切换造成不正确的整体符号翻转。

四元数 $q$ 和 $-q$ 描述同一个 Stokes 旋转，但对应相反的 Jones 算符。由于冻结指标需要 Jones 传播的相位信息，不能为了方便而随意把每步四元数都强制改成 $w\ge0$。代码保留由指数和有序乘法给出的连续算符符号。

当单步横向场轴的净转角超过 0.2 rad 时，代码回退到固定轨道屏幕中的积分，避免在场方向快速转动时依赖转动基底的数值优势。0.2 是数值选择，不是新的物理尺度。

### 13.6 为什么整步和两个半步的差要除以 15

设四阶方法从同一个初值出发，一整步的局部误差为 $Ch^5+O(h^6)$。两个半步的主误差近似相加：

\[
2C(h/2)^5={Ch^5\over16}.
\]

于是

\[
S_{\rm whole}-S_{\rm half}
\approx Ch^5-{Ch^5\over16}={15Ch^5\over16}.
\]

被接受的两个半步结果相对真解的误差约为 $Ch^5/16$，因此估计值为

\[
\boxed{\epsilon\approx{\|S_{\rm half}-S_{\rm whole}\|\over15}.}
\]

代码接受两个半步的结果，不做额外 Richardson 外推，也不通过强制归一化掩盖误差。

局部误差随 h 的五次方变化。若当前误差为 $\epsilon$，目标为 tol，则令

\[
\epsilon(h_{\rm new}/h)^5\approx\mathrm{tol}
\]

可得步长因子 $(\mathrm{tol}/\epsilon)^{1/5}$。程序再乘安全因子 0.9，并把因子限制在 0.2 到 3 之间。

<a id="section-14"></a>

## 14. 仅保存 Stokes 时，如何实现论文的复振幅冻结判据

### 14.1 论文比较的不是 Stokes 差，也不是完整真空载波

FD11 式 (35) 比较慢包络

\[
|\Delta A|=\|A(s+h)-A(s)\|,
\]

并在 $|\Delta A|r/h<10^{-3}$ 时认为冻结。这里的 A 已经去掉了 $e^{ik_0s}$ 的真空载波；若对完整载波电场求差，即使没有任何磁场也会一直快速振荡，那就不是论文的判据。

但是 A 还包含传播矩阵迹部分 $h_0$ 产生的公共相位。Stokes 把这部分相位完全消去了，所以仅比较 $|\Delta S|$ 无法复现式 (35)。

例如固定方向的纯 O 模只获得相位：$A(s+h)=e^{i\lambda_Oh}A(s)$。它的 Stokes 完全不变，但

\[
|\Delta A|=|e^{i\lambda_Oh}-1|
=2|\sin(\lambda_Oh/2)|
\]

一般并不为零。

### 14.2 不需要知道初始公共相位，只需要当前步的完整算符

设当前步的传播算符为

\[
T=e^{i\Phi}G,
\qquad
G=w\mathbf1+i(v\cdot T_{\rm Pauli}),
\]

其中 $T_{\rm Pauli}=(T_Q,T_U,T_V)$ 是第 12.3 节定义的矩阵三元组，避免与总传播算符 T 混淆。

因为未知的初始公共相位 $e^{i\chi_0}$ 会同时乘到当前步的两个端点上，

\[
\|T(e^{i\chi_0}A)-e^{i\chi_0}A\|
=\|e^{i\chi_0}(T-\mathbf1)A\|
=\|(T-\mathbf1)A\|.
\]

所以不用保存“从光子发射到现在累计了多少公共相位”。只要用规定的 H 计算当前步的 $\Phi=\int h_0ds$，就能恢复这一差值的模。

这个结论只对任意**常数初相位**成立。若任意改变沿路的包络规范 $A\to e^{i\chi(s)}A$，冻结指标会改变。因此必须沿用 FD11 的包络约定并保留该约定下的 $h_0$，不能擅自把传播矩阵改成无迹后仍宣称实现了相同的振幅冻结指标。

### 14.3 从范数平方逐步消去复振幅

对归一化纯态 $A^\dagger A=1$，

\[
\begin{aligned}
\|(T-\mathbf1)A\|^2
&=A^\dagger(T^\dagger-\mathbf1)(T-\mathbf1)A\\
&=A^\dagger(T^\dagger T-T^\dagger-T+\mathbf1)A.
\end{aligned}
\]

传播无吸收，$T^\dagger T=\mathbf1$，所以

\[
\|(T-\mathbf1)A\|^2
=2-A^\dagger T^\dagger A-A^\dagger TA
=2\left[1-\operatorname{Re}(A^\dagger TA)\right].
\]

由 $A^\dagger T_iA=S_i$，有

\[
A^\dagger GA=w+i(v_QQ+v_UU+v_VV)=w+i\,v\cdot S.
\]

记 $d_S=v\cdot S$。乘上公共相位并展开：

\[
\begin{aligned}
e^{i\Phi}(w+id_S)
&=(\cos\Phi+i\sin\Phi)(w+id_S)\\
&=(w\cos\Phi-d_S\sin\Phi)
+i(w\sin\Phi+d_S\cos\Phi).
\end{aligned}
\]

因此

\[
\boxed{|\Delta A|^2
=2\left[1-w\cos\Phi+(v\cdot S)\sin\Phi\right].}
\]

右边只有当前光子的三个实 Stokes，以及本步临时生成的传播算符。正式程序没有引入第二份需要持久演化的复电场状态。

### 14.4 小相位下的稳定计算

冻结附近 $w\approx1,\cos\Phi\approx1$，直接计算 $1-w\cos\Phi$ 可能有严重相消。先分解为

\[
1-w\cos\Phi=(1-w)+w(1-\cos\Phi).
\]

由四元数归一化，

\[
(1-w)(1+w)=1-w^2=|v|^2,
\]

所以在 $w>0$ 时可稳定地写成 $1-w=|v|^2/(1+w)$。再用

\[
1-\cos\Phi=2\sin^2(\Phi/2),
\]

得到实际计算式

\[
\boxed{|\Delta A|^2
=2\left[{|v|^2\over1+w}+2w\sin^2(\Phi/2)
+(v\cdot S)\sin\Phi\right]\quad(w>0).}
\]

当 $w\le0$ 时，$1-w$ 本身不接近零，代码直接计算它，以避免 $1+w$ 接近零。理论值非负；最后的 `max(0,d2)` 只防止舍入造成很小的负数。

### 14.5 微小步长极限：独立检查公共相位系数

当 h 很小时，

\[
A(s+h)-A(s)=ihHA+O(h^2),
\]

所以

\[
\lim_{h\to0}{|\Delta A|\over h}
=\sqrt{A^\dagger H^2A}.
\]

写 $H=h_0\mathbf1+dT_Q+gT_U$。由于 $T_Q^2=T_U^2=\mathbf1$ 且反对易子为零，

\[
H^2=(h_0^2+d^2+g^2)\mathbf1
+2h_0(dT_Q+gT_U).
\]

取期望值：

\[
A^\dagger H^2A=h_0^2+d^2+g^2+2h_0(dQ+gU).
\]

代入 $h_0=11\kappa/6$、$d^2+g^2=\kappa^2/4$，以及

\[
dQ+gU={\kappa\over2}(Q\cos\Theta+U\sin\Theta)
={\kappa\over2}Q_B,
\]

得到

\[
\begin{aligned}
A^\dagger H^2A
&=\kappa^2\left({121\over36}+{1\over4}+{11\over6}Q_B\right)\\
&=\kappa^2{65+33Q_B\over18}.
\end{aligned}
\]

因此

\[
\boxed{\lim_{h\to0}{r|\Delta A|\over h}
=\kappa r\sqrt{{65+33Q_B\over18}}.}
\]

O 模 $Q_B=1$ 时为 $7\kappa r/3$，E 模 $Q_B=-1$ 时为 $4\kappa r/3$。这说明即使纯模的 Stokes 恒定，论文的振幅冻结指标依然要等双折射相位足够弱才趋近零。

### 14.6 代码中最终使用的冻结判断

每个接受的步使用该步的合成算符计算

\[
\eta_{\rm freeze,num}={r_{\rm end}\over h}|\Delta A|.
\]

当它小于默认阈值 $10^{-3}$，并且光子向外传播、终点满足 $\kappa r<1$ 时，停止继续求解双折射。

其中需要区分三个层次：

1. $|\Delta A|r/h<10^{-3}$ 是 FD11 的冻结指标和阈值。
2. 取步的终点半径 $r_{\rm end}$ 是当前有限步长实现的约定；对向外传播，它比取起点半径稍保守。
3. 向外传播及 $\kappa r<1$ 是程序附加的数值保护。强相位区可能恰好绕过整数个 $2\pi$，使某一步的端点振幅差很小；这并不表示后续已经冻结，因此需要排除此种伪触发。

停止的是双折射导致的 Stokes 演化。之后仍沿原 GR 光路传播，在平行输运屏幕中保存现有数值；如果再次散射，则重新抽取模式并重新判断启动条件。

<a id="section-15"></a>

## 15. 径向纯模光子的解析解与快速路径

### 15.1 径向光路上为什么磁场投影方向不变

原轨道方程在 $\alpha=0$ 时给出

\[
\psi'=0,\qquad \alpha'=0,\qquad r'=L(r).
\]

因此光子一直沿相同的径向单位向量传播，轨道屏幕也不发生屏幕内转动。

原磁场是自相似场，可以写成

\[
\boldsymbol B(r,\vartheta)
=\left({R_*\over r}\right)^{2+p}\boldsymbol F(\vartheta),
\]

其中 $\boldsymbol F$ 包含固定磁余纬处的三个角向分量。径向光路上的磁余纬不变，所以磁场三个分量只共同乘上一个随 r 变化的实尺度。其横向方向以及 $\Theta=2\phi_B$ 都不变。

### 15.2 推导代码中的径向 $\kappa$ 缩放

比较 $r$ 与参考位置 $r_0$：

\[
B_\perp^2(r)=B_\perp^2(r_0)
\left({r_0\over r}\right)^{2(2+p)}.
\]

又因 $E_\infty$ 沿自由传播段不变，

\[
{E_{\rm loc}(r)\over E_{\rm loc}(r_0)}
={L(r_0)\over L(r)}.
\]

而 $\kappa\propto E_{\rm loc}B_\perp^2$，所以

\[
\boxed{\kappa(r)=\kappa(r_0)
\left({r_0\over r}\right)^{2(2+p)}{L(r_0)\over L(r)}.}
\]

有了参考位置的横向场方向和 $\kappa(r_0)$，后续不必反复查磁场角度表。这就是 `radial_coeff` 缓存的数学依据。

### 15.3 “径向纯模恒定”不是“所有径向偏振都恒定”

径向光路上 $\Theta'=0$，转动基底方程退化为

\[
Q_B'=0,\qquad U_B'=-\kappa V,\qquad V'=\kappa U_B.
\]

对从 O/E 模开始的光子，$U_B=V=0$，右边全部为零，所以它们始终保持 $Q_B=\pm1$。变回轨道屏幕后

\[
S=\pm(\cos\Theta,\sin\Theta,0)
\]

恒定，这是一个精确的 Stokes 解析解。

但是若人为给径向光子一个非本征初态，使 $U_B$ 或 V 不为零，它仍会在 $(U_B,V)$ 平面上转动。当前解析快速路径适用于程序物理流程中由发射或散射产生的径向**纯模**，并不是任意径向混合偏振的通用替代。

### 15.4 纯模的相位与冻结指标

由本征值

\[
\lambda_O=h_0+{\kappa\over2}
=\left({11\over6}+{3\over6}\right)\kappa={7\over3}\kappa,
\]

\[
\lambda_E=h_0-{\kappa\over2}
=\left({11\over6}-{3\over6}\right)\kappa={4\over3}\kappa.
\]

令 $c_O=7/3,c_E=4/3$，当前步的振幅相位增量为

\[
\Delta\varphi_{O/E}=c_{O/E}\int_s^{s+h}\kappa(s')ds'.
\]

纯模振幅满足 $A(s+h)=e^{i\Delta\varphi}A(s)$。因此

\[
\begin{aligned}
|\Delta A|^2
&=|e^{i\Delta\varphi}-1|^2\,|A(s)|^2\\
&=(\cos\Delta\varphi-1)^2+\sin^2\Delta\varphi\\
&=2-2\cos\Delta\varphi
=4\sin^2(\Delta\varphi/2),
\end{aligned}
\]

故

\[
\boxed{|\Delta A|=2\left|\sin\left(
{c_{O/E}\over2}\int_s^{s+h}\kappa(s')ds'\right)\right|.}
\]

代码仍按式 (34) 决定何时进入偏振演化阶段；进入后用精确 Stokes 解和上面的相位积分检测冻结。省去的是不必要的向量 ODE 数值求解，没有省掉论文的冻结判据。

相位积分用两个半区间各自的两点 Gauss 公式。每个半区间的长度为 $h/2$，因此四个节点在全步中的相对位置为

\[
{1\over4}\pm{\gamma_G\over2},\qquad
{3\over4}\pm{\gamma_G\over2},
\]

四个权重均为 $h/4$，对应代码中的四点循环。

<a id="section-16"></a>

## 16. 系数插值、求导与稠密偏振查询

**实现更新：** `fd11.cpp` 已按要求删除系数多项式插值。以下涉及八节点多项式、差商、
Horner 求值和插值误差检查的内容仅保留为历史推导，不再描述当前代码。
当前系数统一由 `pol_coeff` 直接计算，方向导数使用第 7 节的局部角差有限差分；
稠密偏振查询仍保留。删除范围与验证结果见 [fd11_remove_polynomial.md](fd11_remove_polynomial.md)。

### 16.1 为什么插值 $\Omega_Q,\Omega_U$，而不直接插值角度

角度 $\phi_B$ 或 $\Theta$ 有分支和周期性。直接插值它们可能把等价的 $\pi$ 与 $-\pi$ 当成相差很大的两个数。

代码实际插值

\[
v(s)=(\Omega_Q,\Omega_U,0)
=(\kappa\cos\Theta,\kappa\sin\Theta,0),
\]

再用

\[
\kappa=\sqrt{v_x^2+v_y^2},\qquad
\cos\Theta={v_x\over\kappa},\qquad
\sin\Theta={v_y\over\kappa}
\]

恢复系数。这些笛卡尔分量通常比角度本身更适合插值。

### 16.2 八个等距节点与 Newton 差商

设当前轨道段为 $[s_L,s_R]$，宽度 $H_s=s_R-s_L$。定义无量纲变量

\[
t={7(s-s_L)\over H_s}.
\]

八个节点 $s_j=s_L+jH_s/7$ 对应整数 $t_j=j$，$j=0,\ldots,7$。先在这些位置计算 $v_j=v(s_j)$。

零阶差商就是节点值：$v[t_j]=v_j$。高阶差商递推为

\[
v[t_{i-m},\ldots,t_i]
={v[t_{i-m+1},\ldots,t_i]-v[t_{i-m},\ldots,t_{i-1}]
\over t_i-t_{i-m}}.
\]

由于节点间隔在 t 坐标中为 1，分母 $t_i-t_{i-m}=m$。这解释了程序内层循环为什么除以 `order`，而不是再除以一个物理长度。

最后 Newton 多项式为

\[
P(t)=a_0+a_1t+a_2t(t-1)+\cdots
+a_7t(t-1)\cdots(t-6),
\]

其中 $a_j=v[t_0,\ldots,t_j]$。向量三个分量分别使用同一套差商操作。

### 16.3 Horner 递推与导数

从最高阶开始令 $P_7=a_7$，然后递推

\[
P_i=a_i+(t-i)P_{i+1},\qquad i=6,\ldots,0.
\]

这就是代码的 Horner 求值。对每一步求导：

\[
{dP_i\over dt}=P_{i+1}+(t-i){dP_{i+1}\over dt}.
\]

所以可在同一个向后循环中同时更新向量值与导数。最后用链式法则

\[
{dP\over ds}={dP\over dt}{dt\over ds}={7\over H_s}{dP\over dt}.
\]

把这个导数代入第 7.3 节的 $\Theta'$ 公式，就得到代码中的前因子 `7/poly_h`。

### 16.4 插值误差检查为何还要乘上轨道段长度

额外采样点上，代码比较插值与直接磁场计算的 $\Omega$。除相对误差限制外，还检查

\[
|\delta\Omega|H_s\lesssim0.01\,\mathrm{tol}.
\]

其动机是：相位来自生成元沿路的积分。若整段都有统一误差上界 $\epsilon_\Omega$，则

\[
\left|\int\delta\kappa\,ds\right|
\le H_s\sup|\delta\kappa|
\le H_s\sup|\delta\Omega|.
\]

强双折射区即使相对系数误差很小，也可能积累明显相位误差，因此只检查相对误差不够。

需要强调：程序只在额外的几个采样点验证，**没有证明**这些点给出了全区间误差的严格上界。该条件是数值控制策略，最终精度还需要独立参考解和容限收敛测试。检查失败时，代码回到直接磁场求值。

### 16.5 光学深度积分怎样查询任意位置的 Stokes

每个接受的偏振半步保存其起点位置、终点位置和起点 Stokes。求积器查询半步内部的某点 $s_q$ 时，先找到包含它的分段，然后从该段起点重新构造一个长度为 $s_q-s_{\rm begin}$ 的传播算符，作用在该段起始 Stokes 上。

因此这里的“稠密偏振解”不是把 Q、U、V 分别作普通线性插值，而是在已接受的小区间内重新应用同一传播算法。光学深度积分会在这些段的连接点处分段，避免求积器在容差量级的插值接缝处反复细分。

<a id="section-17"></a>

## 17. 逃逸时的统一屏幕与输出

### 17.1 为什么还要做一次坐标变换

每个光子的轨道平面一般不同，所以它的内部 $e_x=n\times k,e_y=n$ 也不同。直接把不同光子的内部 Q、U 相加，会混合不同的参考方向。

代码选择磁轴方向 $\hat m=(0,0,1)$，把它投影到垂直传播方向的平面：

\[
m_\perp=\hat m-(\hat m\cdot\hat k)\hat k.
\]

直接检查

\[
m_\perp\cdot\hat k=\hat m\cdot\hat k
-(\hat m\cdot\hat k)(\hat k\cdot\hat k)=0.
\]

若 $m_\perp\ne0$，定义输出基

\[
e_{x,\rm out}={m_\perp\over|m_\perp|},\qquad
e_{y,\rm out}=\hat k\times e_{x,\rm out}.
\]

它也是右手横向基。若光子方向恰好沿磁轴，投影为零，此时输出参考轴在几何上不唯一，代码固定回退到原轨道屏幕 x 轴。

### 17.2 用点积计算输出角度

设输出 x 轴相对内部屏幕 x 轴的角度为 $\chi$。由于两个屏幕正交归一，

\[
c_\chi=e_{x,\rm out}\cdot e_x=\cos\chi,
\qquad
s_\chi=e_{x,\rm out}\cdot e_y=\sin\chi.
\]

无需显式求角度，直接构造

\[
\cos2\chi=c_\chi^2-s_\chi^2,\qquad
\sin2\chi=2c_\chi s_\chi.
\]

再使用第 7.1 节已经完整推导的变换：

\[
\boxed{
Q_{\rm out}=Q\cos2\chi+U\sin2\chi,\quad
U_{\rm out}=-Q\sin2\chi+U\cos2\chi,\quad
V_{\rm out}=V.}
\]

输出的五列为 $E_\infty,\mu_k,Q_{\rm out},U_{\rm out},V_{\rm out}$，其中 $\mu_k=\hat k\cdot\hat m$，每个光子的 $I=1$。

### 17.3 逃逸不等于偏振已经冻结

沿用原程序的逃逸边界 $r=10000$ km。代码在此处记录当前偏振：

- 若已经冻结，输出保存的 Stokes 经上述屏幕变换后的数值。
- 若正在积分，输出边界处当前积分值。
- 若仍处于 O/E 标签阶段，就在边界处把该本征模转换成 Stokes 后输出，不为此补做一段此前未启动的 ODE。

因此程序输出的是指定边界处的偏振，不自动等同于无穷远处的极限值。

<a id="section-18"></a>

## 18. 与上述推导有关的原有几何公式补充

以下公式不是上一任务新改动的物理模型，但它们进入 $\kappa,\mu$ 的计算。为使从轨道变量到偏振系数的链条完整，在此补出推导。

### 18.1 原轨道方程中的 $\alpha'$

球对称静态时空中，角动量与能量之比给出守恒的冲量参数

\[
b_{\rm imp}={r\sin\alpha\over L(r)}.
\]

在非径向且可直接取对数的区间，对它求导：

\[
0={r'\over r}+\cot\alpha\,\alpha'-{L_r\over L}r',
\]

其中 $L_r=dL/dr$。整理并代入 $r'=L\cos\alpha$：

\[
\alpha'=-\left({1\over r}-{L_r\over L}\right)
L\cos\alpha\tan\alpha
=-L\sin\alpha\left({1\over r}-{L_r\over L}\right).
\]

由 $L^2=1-r_s/r$ 得

\[
2LL_r={r_s\over r^2},\qquad
{L_r\over L}={r_s\over2r^2L^2}.
\]

代入：

\[
\begin{aligned}
\alpha'
&=-{\sin\alpha\over rL}
\left(L^2-{r_s\over2r}\right)\\
&=\boxed{-{\sin\alpha\over rL}
\left(1-{3r_s\over2r}\right)}.
\end{aligned}
\]

这个结果在转向点和径向极限可按连续性延拓；径向极限直接给出 $\alpha'=0$。与第 2.2 节的两个方程合在一起，就是原 `geodesic()` 的三个分量。

### 18.2 从球坐标磁场分量变到屏幕分量

设位置方向写为

\[
\hat r=(\sin\vartheta\cos\varphi,\sin\vartheta\sin\varphi,\cos\vartheta).
\]

分别对角度求导并归一化，可得

\[
\hat\vartheta=(\cos\vartheta\cos\varphi,\cos\vartheta\sin\varphi,-\sin\vartheta),
\]

\[
\hat\varphi=(-\sin\varphi,\cos\varphi,0).
\]

令位置单位向量的笛卡尔分量为 $(r_x,r_y,r_z)$，并记 $\rho=\sqrt{r_x^2+r_y^2}=\sin\vartheta$。非极点处可写为

\[
\hat\vartheta=\left({r_xr_z\over\rho},{r_yr_z\over\rho},-\rho\right),
\qquad
\hat\varphi=\left(-{r_y\over\rho},{r_x\over\rho},0\right).
\]

磁场表返回当地球坐标正交分量 $(B_r,B_\vartheta,B_\varphi)$，因此

\[
\boldsymbol B=B_r\hat r+B_\vartheta\hat\vartheta+B_\varphi\hat\varphi.
\]

最后取点积

\[
B_x=\boldsymbol B\cdot e_x,\qquad
B_y=\boldsymbol B\cdot e_y,\qquad
\mu={\boldsymbol B\cdot\hat k\over B}.
\]

极点处球坐标方位方向本身不唯一，代码选固定正交方向作为退化约定。

### 18.3 原 `geo()` 中的方向余弦表达式

利用 $\hat k=\cos\alpha\hat r+\sin\alpha(n\times\hat r)$，

\[
\mu=b_r\cos\alpha+
\sin\alpha\left[b_\vartheta(n\times\hat r)\cdot\hat\vartheta
+b_\varphi(n\times\hat r)\cdot\hat\varphi\right].
\]

标量三重积给出

\[
(n\times\hat r)\cdot\hat\vartheta
=n\cdot(\hat r\times\hat\vartheta)=n\cdot\hat\varphi,
\]

\[
(n\times\hat r)\cdot\hat\varphi
=n\cdot(\hat r\times\hat\varphi)=-n\cdot\hat\vartheta.
\]

所以

\[
\boxed{\mu=b_r\cos\alpha+
\sin\alpha\left[b_\vartheta(n\cdot\hat\varphi)
-b_\varphi(n\cdot\hat\vartheta)\right],}
\]

与 `geo()` 和 `perform_scattering()` 中的分量形式一致。

<a id="section-19"></a>

## 19. 独立复振幅验证方程与实现对应表

### 19.1 参考测试中的四个实数方程

正式程序只存 Stokes；为了独立验证，测试程序把 Jones 振幅写成

\[
A_x=x_R+ix_I,\qquad A_y=y_R+iy_I.
\]

把它们代入 $A_x'=i(h_0+d)A_x+igA_y$，实部和虚部分别为

\[
x_R'=-(h_0+d)x_I-gy_I,
\qquad x_I'=(h_0+d)x_R+gy_R.
\]

再对第二个 Jones 方程作相同展开：

\[
y_R'=-gx_I-(h_0-d)y_I,
\qquad y_I'=gx_R+(h_0-d)y_R.
\]

参考解用独立的 DOPRI5 积分这四个实数，再通过

\[
Q=x_R^2+x_I^2-y_R^2-y_I^2,
\]

\[
U=2(x_Ry_R+x_Iy_I),\qquad
V=2(x_Iy_R-x_Ry_I)
\]

转成 Stokes，与正式算法对照。这种验证同时保留公共相位，因而也能检查第 14 节的振幅差表达式。

### 19.2 公式与代码的对应

| 本文推导内容 | 主要代码位置/函数 | 来源性质 |
|---|---|---|
| 真空传播矩阵、$\kappa$、双倍角 | `pol_coeff()` | FD11 (18)–(24) 的弱场一阶展开 |
| 三个 Stokes 方程、叉乘生成元 | `fd11::Coeff::omega()` | 由 FD11 (36)–(39) 严格代数变换 |
| 绝热判据启动、O/E 初值 | `adiabatic_event()`、`step_geodesic()`、`evolve_geodesic()`、`start_polarization()`、`eigenmode()` | `doc/adiabatic.md` 的事件判据及本征向量转换 |
| 转动磁场基底 | `pol_step()` | 对同一 Stokes ODE 的坐标变换 |
| 旋转和四元数合成 | `Rotation`、`compose()` | 长度守恒传播算符的实数表示 |
| 四阶指数步 | `fd11::magnus()` | 数值积分选择，不是 FD11 原文算法 |
| 步长减半误差与控制 | `advance_polarization()` | 四阶局部截断误差推导 |
| 公共相位恢复、式 (35) 指标 | `Propagator::amplitude_change()` | 对论文振幅差的 Stokes/算符重写 |
| 径向缓存和解析解 | `pol_coeff()`、`advance_polarization()` 径向分支 | 自相似磁场与径向纯模的解析性质 |
| 式 (33) 偏振重叠 | `fd11::overlap()`、`polarized_rates()` | 复振幅模平方的严格展开 |
| $c_0+c_QQ_B+c_VV$ | `opacity_terms()` | 原共振核与新偏振重叠的系数整理 |
| 快速振荡分部积分 | `fast_optical_depth()` | 利用 Stokes 方程的积分恒等式 |
| 支集、共振切点与事件抽样 | `evolve_geodesic()` | 原输运流程，加偏振分段处理 |
| 散射后重置模式 | `perform_scattering()` | 沿用的共振散射模型 |
| 输出屏幕旋转 | `evolve_geodesic()` 的逃逸分支 | 同一 Stokes 状态的参考轴变换 |
| 四实数 Jones 参考方程 | `src/dev/fd11_test.cpp` | 测试用途，不进入正式光子状态 |

### 19.3 理论、近似和数值检查的边界

本文中的 Stokes 变换、旋转基底变换、散射重叠展开、传播算符冻结公式和分部积分恒等式，都是在给定模型下的代数关系。

弱场真空张量、忽略纵向反馈、延迟到式 (34) 才从纯模开始积分、以及散射后抽取纯 O/E 模，是物理模型或论文采用的近似。特别是式 (34) 规定的 $10^{-3}$ 是一个初始化选择，不能自动解释为每条特殊光路最终 Stokes 都具有 $10^{-3}$ 的物理精度。

轨道采样、有限差分、七次插值、四阶积分容限、有限步长冻结指标和有限逃逸半径，则产生数值或边界误差。这些误差需要参考解、加密或参数敏感性检查；单看 $Q^2+U^2+V^2=1$ 不足以证明偏振方向正确，因为一个错误方向的旋转也能保持长度。

上一任务的参考解、误差、性能和参数敏感性结果统一记录在 [fd11.md](fd11.md)。本文件只补充、展开推导，没有修改 `src/dev/fd11.cpp` 或原输运程序。
