#!/usr/bin/env python3
"""Summarize saved measurements without rerunning or changing the reference."""
import csv
import json
from pathlib import Path
import statistics

OUT=Path(__file__).resolve().parent
summary=json.loads((OUT/'summary.json').read_text())

def pick(method,tolerance=1e-6,parameter=.01,initial=.01):
    return next(s for s in summary if s['method']==method and float(s['tolerance'])==tolerance and float(s['parameter'])==parameter and float(s['initial_step'])==initial)

def rows(s):
    with (OUT/s['file']).open() as f:
        next(f)
        return list(csv.DictReader(f))

def save(name,data):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=data[0].keys());w.writeheader();w.writerows(data)

labels={'baseline':'原 DOPRI5','ft07':'FT07 限制器 + DOPRI5（截断 99.8%）','thermal_cap':'温度 log(x) 限制器 + DOPRI5','spatial_regularized':'空间正则化 + 预置区间','velocity':'显式支撑区间的速度积分'}
primary=[]
for name in labels:
    s=pick(name)
    primary.append(f"| {labels[name]} | {float(s['batch_median_ms']):.3f} | {s['max_relative']:.6g} | {s['p95_relative']:.3g} | {s['above_1e3']} | {s['mean_evaluations']:.1f} |")

convergence=[]
for name in labels:
    for tol in [1e-6,1e-7,1e-8,1e-10]:
        s=pick(name,tol)
        convergence.append(f"| {name} | {tol:.0e} | {float(s['batch_median_ms']):.3f} | {s['max_relative']:.6g} | {s['missed_90pct']} | {s['failures']} |")

with (OUT/'ft07_truncation.csv').open() as f:trunc=list(csv.DictReader(f))
decomposition=[]
for tol in [1e-6,1e-8,1e-10]:
    numeric=rows(pick('ft07',tol))
    for n,t in zip(numeric,trunc):
        decomposition.append(dict(tolerance=tol,**{k:n[k] for k in ['b0','muz','pol','omega_inf']},
            tau=n['tau'],truncated_reference=t['truncated'],
            relative_to_truncated=abs(float(n['tau'])-float(t['truncated']))/float(t['truncated']),
            truncation_fraction=t['omitted_fraction']))
save('ft07_error_decomposition.csv',decomposition)
decomp_max={tol:max(d['relative_to_truncated'] for d in decomposition if d['tolerance']==tol) for tol in [1e-6,1e-8,1e-10]}

best=rows(pick('spatial_regularized',1e-10))
discrepancies=sorted([r for r in best if float(r['rel_error'])>1e-8],key=lambda r:float(r['rel_error']),reverse=True)
save('reference_discrepancies.csv',discrepancies)
worst=discrepancies[0]

stability=[]
for name in ['spatial_regularized','velocity']:
    tight=rows(pick(name,1e-10))
    for tol in [1e-6,1e-7,1e-8]:
        coarse=rows(pick(name,tol))
        stability.append(dict(method=name,tolerance=tol,relative_change_to_1e10=max(abs(float(a['tau'])-float(b['tau']))/abs(float(b['tau'])) for a,b in zip(coarse,tight))))
save('quadrature_convergence.csv',stability)

v=dict(line.split('=',1) for line in (OUT/'validation_summary.txt').read_text().splitlines())
thermal=pick('thermal_cap',1e-8);spatial=pick('spatial_regularized',1e-8)
env=json.loads((OUT/'environment.json').read_text())

report=f'''# 共振层漏积分：结果与建议

已完成五种方法的独立 C++ benchmark 和 31 组参数配置；每组均与原始 900 行参考表逐行比较。**当前向外径向程序推荐使用空间正则化与预置区间方法，已接入 `src/main.cpp`。** 原轨道方程移至 `src/photon_evolution.hpp`，原始算法可由 `benchmark/baseline.cpp` 重现。一般非径向传播、南半球双共振分支以及多次进出共振层不在本次验证范围内。

## 主要结果

以下均为 nominal tolerance=10⁻⁶，FT07/温度限制器参数 0.01、初始 proper-path 步长 0.01 km。DOPRI 使用 `atol=rtol=tol`；求积使用 `atol=0.1 tol, rtol=tol`，这不是完全相同的误差范数，因此下面另列实测精度相近的比较。时间是预热后 **7 批、每批全部 900 例** 的中位数，不是单例耗时，单位 ms。运行平台 `{env['machine']}`，GCC 16.2.0，`-O3 -std=c++23`，无 fast-math。具体环境与输入 SHA256 见 [environment.json](environment.json)。

| 方法 | 900 例耗时/ms | 最大相对误差 | P95 相对误差 | 误差 >0.1% 的例数 | 平均函数评估数 |
|---|---:|---:|---:|---:|---:|
{chr(10).join(primary)}

上述五组均无程序异常。函数评估的含义因方法不同而不同，不能当作等价成本：DOPRI 是四维 RHS，求积是 integrand，后者不包含根求解工作；实际时间包含每例几何准备和根求解。

在实测精度相近时，温度限制器 `tol=10⁻⁸` 的最大误差 {thermal['max_relative']:.6g}、耗时 {float(thermal['batch_median_ms']):.3f} ms；空间正则化 `tol=10⁻⁸` 的最大误差 {spatial['max_relative']:.6g}、耗时 {float(spatial['batch_median_ms']):.3f} ms，约快 **{float(thermal['batch_median_ms'])/float(spatial['batch_median_ms']):.1f} 倍**。这部分优势也来自利用径向几何的不变量，以及无需继续积分共振层外至 10000 km 的零 opacity 区域，不能全部归因于步长策略。

## 漏层与 FT07 对照

原程序 `tol=10⁻⁶` 下有 27 例低估超过 90%。例如 β₀=−0.1、μz=0、E 模、ω∞=0.01：参考 τ=10.090418511430368，原程序 τ≈1.2832×10⁻²⁴。所有 stage 都未采到层内时，嵌入误差接近零，局部误差控制不会自动发现遗漏。

FT07 是最先实现和运行的修复对照。按论文 §3.4 保留中央 99.8% 的粒子，按式 (39) 限制动量相对步长、按式 (31) 转换为空间步长、按式 (38) 限制根合并附近的步长；光深积分仍用原 DOPRI5，物理归一化和 GR 修正沿用 note。具体适配和经验规则见 [ALGORITHMS.md](ALGORITHMS.md)。它不是论文完整 Monte Carlo 程序的复刻。

对完全相同的截断速度区间独立积分后，发现**只截掉 0.2% 粒子即可损失最多 {float(v['max_FT07_truncated_fraction'])*100:.6f}% 光深**，因为散射权重依赖速度。FT07 数值结果相对该截断积分的最大误差分别为：

| DOPRI tolerance | 相对截断区间真值的最大误差 |
|---|---:|
| 10⁻⁶ | {decomp_max[1e-6]:.6g} |
| 10⁻⁸ | {decomp_max[1e-8]:.6g} |
| 10⁻¹⁰ | {decomp_max[1e-10]:.6g} |

因此不能通过收紧 DOPRI 容差消除 FT07 在完整分布基准上的约 1.23% 误差。逐例分解见 [ft07_error_decomposition.csv](ft07_error_decomposition.csv)。

## 容差、步长和失败案例

| 方法 | tolerance | 耗时/ms | 成功例的最大相对误差 | 低估 >90% | 异常例数 |
|---|---:|---:|---:|---:|---:|
{chr(10).join(convergence)}

在 `tol=10⁻¹⁰`，原 DOPRI 和温度限制器均有 22 个赤道 E 模算例因 `stepsize underflow` 失败。这是空间 opacity 的可积奇点造成的另一种困难。失败例保留在 CSV 中，**不算入“成功例最大误差”**；不能据剩余算例的精度声称该配置整体成功。

初始步长从 0.01 改为 0.001、0.1、1 km 后，原程序严重漏层例数分别为 28、26、21，说明结果依赖 stage 采样相位；温度限制器这三组均无严重漏层，最大相对误差不超过 7.59×10⁻⁴。FT07 分别测试了 0.005/0.01/0.02 的速度限制系数；温度限制器测试了 0.002/0.01/0.05/0.1。较细步长不保证表观误差严格单调，尤其在截断边界和可积奇点附近；全部配置见 [summary.csv](summary.csv)。

## 推荐方法与精度边界

首先定位共振支撑边缘 `x(r_end)=1`，用 `r=r_end(1−u²)` 消除赤道平方根奇点，并在分布核心对应的空间位置预先分段，再执行自适应 Gauss–Legendre 求积。这样同时解决“没有采到层”和“奇点迫使步长下溢”。完整速度分布保留，`integrate_to(r)` 可用于求累计光深与散射半径。

空间正则化和独立速度积分在 900 例上的最大相对差为 {float(v['max_spatial_velocity_relative']):.3g}；但它们相对用户参考表的最大误差仍为 **{float(worst['rel_error']):.8g}**。有 {len(discrepancies)} 例对表相对误差超过 10⁻⁸，主要集中在高能端。这些差异已全部计入指标，**没有改表、删除算例或用新积分结果替代参考值**。

最明显的一例：β₀=−0.9、μz=0.4、O 模、ω∞=100，表中 τ={float(worst['reference']):.17g}，空间积分 τ={float(worst['tau']):.17g}。收紧本次求积容差没有消除该差异；本报告将它视为对给定参考的残余偏差，不据交叉一致性断言参考表错误，也不声称对表达到 10⁻¹⁰ 精度。详见 [reference_discrepancies.csv](reference_discrepancies.csv) 和 [quadrature_convergence.csv](quadrature_convergence.csv)。

## 验证与复现

除完整表比较外，独立速度积分用于交叉检查：

- 总光深：900 例全部通过。
- 累计光深：{v['cumulative_cases']} 例，最大相对差 {float(v['max_cumulative_relative']):.3g}。
- 固定 U∈{{0.01,0.1,0.5,0.9}} 的散射/逃逸：{v['event_cases']} 例，其中 {v['escapes']} 例逃逸，判断全部一致；发生散射时最大半径相对差 {float(v['max_event_radius_relative']):.3g}。
- 额外冷/高温参数 β₀∈{{−0.03,−0.05,−0.99}}：{v['stress_cases']} 例，最大相对差 {float(v['max_stress_relative']):.3g}。
- 新增 `do_step(max_step)`：零 RHS、控制器自动放大后、正反向积分和非法上限检查通过。
- `src/main.cpp` 编译运行完成，900 行参数顺序与参考一致；输出保存在 [main_results.txt](main_results.txt)。编译启用 `-Wall -Wextra -pedantic`，最终无编译警告。

这些累计量/散射半径验证使用独立公式的 C++ 速度积分，原表本身只有总光深，因此不能说它们直接由原表验证。

完整重跑：

```sh
python3 output/run_benchmarks.py --repeats 7
python3 output/write_report.py
```

单独运行某个方法（从仓库根目录）：

```sh
g++-16 -O3 -std=c++23 benchmark/spatial_regularized.cpp -o output/spatial_regularized
output/spatial_regularized output/my_result.csv 1e-8 .01 7 .01
```

参数顺序为：输出 CSV、容差、限制器系数、重复次数、初始步长。后两类径向求积算法忽略限制器系数和初始步长；统一命令行只是方便比较。`Makefile` 原有路径指向不存在的 `ode/main.cpp`，本次没有改动根目录文件；请使用上述已验证命令或 output/ 中的脚本。

所有新增方法分别在 `benchmark/*.cpp`，生产算法在 `src/`，结果、辅助测试文件、二进制和说明均在 `output/`。完整算法推导、FT07 的实现细节及计时范围见 [ALGORITHMS.md](ALGORITHMS.md)。
'''
(OUT/'REPORT.md').write_text(report)
print('Wrote output/REPORT.md and diagnostic CSV files.')
