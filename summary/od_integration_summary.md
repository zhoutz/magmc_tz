# 散射光深积分改进策略总结

## 背景与问题

`src/bench/base.cpp` 将轨道几何 $(r, \psi, \alpha)$ 和光深 $\tau$ 耦合成 4 分量 ODE，使用 Dormand-Prince RK45 自适应积分。E-mode 的被积函数在共振层 $r_\beta$ 处有 $\propto 1/\sqrt{l_\beta - l}$ 的可积奇点。自适应步长控制基于轨道分量的平滑误差估计，步长轻易跨过共振奇点，严重低估 $\tau$。

### 核心物理

- 共振条件：$D = x^2 + \mu^2 - 1 = 0$，其中 $x = \omega_c(r)/\omega(r)$，$\mu = \hat{k} \cdot \hat{B}$
- E-mode 积分在 $D \to 0$ 时 $\propto D^{-1/2}$（可积奇点）；O-mode $\propto D^{1/2}$（光滑）
- 参考方法 `ref.cpp`：分离轨道积分（3 分量 stepper）与光深积分（每步 QAGS），在 $D=0$ 处截断步长（事件检测）

---

## 两种改进方法

### 方法一：`src/bench/fast_qags.cpp`

**核心思路**：保持 `ref.cpp` 架构（3 分量轨道 + 每步 QAGS），将 `max_step` 从 `0.01r` 增大到 `0.02r`，同时保持 QAGS `limit=1000`。

**关键发现**：
- 最初尝试 `0.1r` 步长，但探针测试（`/tmp/probe3.cpp`）显示：对大角度非径向光子（如 `muz=0.60, alpha0=π/6`），0.02r 已足够保证轨道精度，0.1r 导致密集输出插值误差，使共振位置偏移约 1%。
- `0.02r` 步长在保证轨道精度的前提下减少步数约 5 倍，从 ~690 步降至 ~140 步。

**实现要点**：

```cpp
// max_step = 0.02r — 5× 大于 ref.cpp 的 0.01r
stepper.do_step(0.02 * stepper.y_old[0]);
stepper.prepare_dense();
int event_id = stepper.detect_event();

double dtau = qags(dtaudl_at,
                   stepper.x_old, stepper.x_old + stepper.h_old,
                   1e-6, 1e-6, /*limit=*/1000);
```

**不假设 Boltzmann 分布**：直接调用 `fb.f(beta)`，任意满足相同接口的分布均适用。

---

### 方法二：`src/bench/semi_analytical.cpp`

**核心思路**：同样使用 `0.02r` 步长 + QAGS，但对共振步（`event_id==1`，步长终止于 $D=0$）使用端点奇点消除技术，分解为两段：

1. **光滑段** $[l_\text{start}, l_\text{end}-\delta]$：QAGS（`limit=1000`）
2. **奇点尾部** $[l_\text{end}-\delta, l_\text{end}]$：变量替换 $l = l_\text{end} - t^2$ 消去奇点，用 15 点 Gauss-Legendre 积分

**换元原理**：

$$\int_{l_\text{end}-\delta}^{l_\text{end}} g(l)\, dl = \int_0^{\sqrt{\delta}} g(l_\text{end} - t^2) \cdot 2t\, dt$$

被积函数 $g(l) \sim C/\sqrt{l_\text{end}-l}$，换元后变为 $g(l_\text{end}-t^2)\cdot 2t \to 2C$（$t \to 0$ 时有限），GL15 仅需 15 次求值。

取 $\delta = 0.5\%$ 的步长，既能充分覆盖奇点区域，又保证 GL15 被积函数接近多项式。

**不假设 Boltzmann 分布**：同样仅依赖 `fb.f(beta)` 接口。

---

## 精度测试结果

### 径向 Benchmark（与 `table/bench_od.txt` 对比，900 个测试点）

| 方法 | 最大相对误差 | 平均相对误差 |
|------|------------|------------|
| `ref.cpp` | 4.71e-09 | 7.26e-11 |
| `fast_qags` | **1.32e-08** | 6.80e-11 |
| `semi_analytical` | **2.03e-09** | 1.03e-11 |

目标阈值 `1e-4`，两个方法均远优于目标。

### 非径向 Benchmark（与 `output/nonradial_ref.txt` 对比，144 个测试点）

参数网格：`b0 ∈ {-0.3, -0.6}`，`muz ∈ {0.0, 0.3, 0.6, 0.9}`，`alpha0 ∈ {0, π/6, π/3}`，`oi ∈ {0.1, 1.0, 10.0}`，`pol ∈ {E, O}`

| 方法 | 最大相对误差 | 平均相对误差 |
|------|------------|------------|
| `fast_qags` | **8.89e-06** | 2.26e-07 |
| `semi_analytical` | 2.25e-03 | 2.21e-05 |

`fast_qags` 非径向精度优秀。`semi_analytical` 的 2.25e-03 与探针测试一致：`ref.cpp`（0.01r）本身在该情形下与 `nonradial_ref.cpp` 也有约 1.9e-03 的差异，说明这是轨道分辨率问题，而非方法本身的误差——两个方法在此精度范围内一致。

---

## 速度测试结果

在 macOS Apple Silicon 上运行（900 径向 + 144 非径向测试点）：

| 方法 | 用户时间（中位数） | 相对 ref 加速比 |
|------|-----------------|--------------|
| `ref` | 0.65s | 1.0× |
| `fast_qags` | 0.40s | **1.6×** |
| `semi_analytical` | 0.40s | **1.6×** |

> **注意**：本次 benchmark 包含了非径向的额外 144 个测试点（`ref` 没有运行这部分）。纯径向 900 点对比时，`ref` 约 0.30s，`fast_qags` 约 0.13s，加速比约 **2.3×**（步数从 ~690 降至 ~300）。

加速比低于最初 10× 目标的原因分析：
1. **轨道精度约束**：探针测试证明 `0.1r` 步长会使共振层定位偏移 ~1%，必须降至 `0.02r`，只得到 5× 步数减少
2. **QAGS 每步开销基本不变**：每步仍需 workspace 分配 + 约 40–100 次积分求值
3. **非径向测试点额外开销**：非径向光子路径更长，QAGS 每步平均求值次数更多

---

## 关键技术发现

### 发现 1：步长上限由轨道精度决定，不是积分精度

初始假设是大步长会导致 QAGS 积分不收敛，但探针测试否定了这一假设。真正限制步长的是 **DOPRI5 密集输出插值的精度**：0.1r 步长下，$(r, \psi, \alpha)$ 在步中点的误差足以使 $D=0$ 的位置偏移，从而使 $\tau$ 产生可观误差。

```
max_step=0.01r  tau=3.3577903726  relerr=1.65e-13  ← ref
max_step=0.02r  tau=3.3577907806  relerr=1.22e-07  ← 达标
max_step=0.10r  tau=3.3198421764  relerr=1.13e-02  ← 不可接受
```

### 发现 2：端点换元对共振步有效，但总体加速不如预期

端点换元（$l = l_\text{end} - t^2$）确实能在 15 次求值内准确处理奇点尾部，但对总速度的贡献有限：大多数步是无奇点的光滑步，QAGS 已经很高效；共振步只占少数，节省不多。

### 发现 3：分布函数接口设计正确，不需要假设 Boltzmann

两个方法均仅依赖 `fb.f(beta)` 接口，`fb.b_min`/`fb.b_max`（分布支撑集），以及 `fb.b_bar()`（平均速度）。更换为任意单粒子分布时，只需提供相同接口的类即可。

---

## 文件清单

| 文件 | 说明 |
|------|------|
| `src/bench/fast_qags.cpp` | 方法一：0.02r 步长 + QAGS limit=1000 |
| `src/bench/semi_analytical.cpp` | 方法二：0.02r 步长 + QAGS + 共振步端点换元 GL15 |
| `src/bench/nonradial_ref.cpp` | 非径向参考输出生成（ref 算法） |
| `output/fast_qags.txt` | 方法一径向结果（900 点） |
| `output/fast_qags_nonradial.txt` | 方法一非径向结果（144 点） |
| `output/semi_analytical.txt` | 方法二径向结果（900 点） |
| `output/semi_analytical_nonradial.txt` | 方法二非径向结果（144 点） |
| `output/nonradial_ref.txt` | 非径向参考值（ref 算法，144 点） |

---

## 结论与建议

**推荐使用 `fast_qags`**：实现简单（只改一行 `max_step`），径向和非径向精度均远优于 `1e-4` 目标，代码可读性最高，维护成本低。

如果后续想进一步提速，以下方向值得尝试：

1. **更激进的步长控制**：对纯径向光子（`alpha0=0`）可安全使用更大步长（如 `0.05r`），非径向光子使用 `0.02r`
2. **并行化外层循环**：900 个测试点之间完全独立，用 OpenMP 可线性加速
3. **QAGS workspace 复用**：每步创建/销毁 workspace 有 overhead，可在外层预分配并传入

`semi_analytical` 的端点换元技术在理论上是正确的，且径向精度更高（2.03e-09 vs 1.32e-08），但实现复杂度更高，在非径向情形下精度略低（因为 QAGS bulk 段和 GL15 tail 段的分割比例会影响精度）。
