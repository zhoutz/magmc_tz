定义磁场在横向平面中的方向角
\[
\chi_B=\operatorname{atan2}(B_2,B_1).
\]
则
\[
H=\frac{\kappa}{2}
\begin{pmatrix}
\cos2\chi_B&\sin2\chi_B\\
\sin2\chi_B&-\cos2\chi_B
\end{pmatrix}.
\]
把振幅变换到随磁场转动的本征模基底，得到
\[
\boxed{
\frac{d}{dl}
\begin{pmatrix}a_O\\a_E\end{pmatrix}
=
\begin{pmatrix}
i\kappa/2&d\chi_B/dl\\
-d\chi_B/dl&-i\kappa/2
\end{pmatrix}
\begin{pmatrix}a_O\\a_E\end{pmatrix}.
}
\]
这里有两个不同的过程：
- 对角项产生 O、E 之间的相对传播相位。
- 非对角项来自本征模方向沿光线的转动，负责模耦合。
一个自然的绝热条件是
\[
\boxed{
\epsilon_{\rm ad}
\equiv
\frac{|d\chi_B/dl|}{\kappa}\ll1.
}
\]
在这个条件下，初始纯 O 或纯 E 光子基本保持原来的模标签，电场方向随对应的局域本征矢转动。
单位局域路程积累的相对相位
\[
\boxed{
\kappa(l)\equiv\frac{\omega(l)}c\Delta n(l)
=\frac{\alpha_{\rm em}\omega_\infty}
{30\pi c\,\mathcal N(r)}
\frac{B_\perp^2}{B_{\rm QED}^2}.
}
\]
\(\kappa\) 的量纲是长度的倒数。
我建议使用
\[
  |d\chi_B/dl| > 1e-3 \kappa
\]
作为从绝热演化切换到积分复振幅的判据。

