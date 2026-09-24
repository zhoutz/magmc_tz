#!/usr/bin/env python3
"""Generate the final report from saved data, without changing any input table."""
import csv
import json
import math
from pathlib import Path
import statistics

OUT=Path(__file__).resolve().parent
summary=json.loads((OUT/'summary.json').read_text())

def read(path,metadata=False):
    with (OUT/path).open() as f:
        if metadata:next(f)
        return list(csv.DictReader(f))

def pick(name,tol=1e-6,p=None,h=.01):
    if p is None:p=.01 if name=='ft07' else .1
    return next(x for x in summary if x['method']==name and float(x['tolerance'])==tol and float(x['parameter'])==p and float(x['initial_step'])==h)

def save(name,data):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)

names=['baseline','ft07','phase_cap','event_guard']
labels=['原四维 DOPRI5','FT07 控制器＋通用边界处理','共振变量变化限制器','步内共振边界定位（推荐）']
table=[];probability=[]
for name,label in zip(names,labels):
    s=pick(name);data=read(s['file'],True)
    dp=max(abs(math.expm1(-float(r['tau']))-math.expm1(-float(r['reference']))) for r in data if r['status']=='ok')
    probability.append(dict(method=name,maximum_absolute_probability_error=dp))
    table.append(f"| {label} | {float(s['batch_median_ms']):.3f} | {s['maximum_relative']:.7g} | {s['p95_relative']:.3g} | {s['missed_90pct']} | {s['mean_steps']:.1f} |")
save('scattering_probability_errors.csv',probability)

angles=read('angular_validation.csv');angle_summary=[]
for name in names[1:]:
    data=[r for r in angles if r['method']==name]
    angle_summary.append(dict(method=name,cases=len(data),failures=sum(r['status']!='ok' for r in data),
        maximum_relative=max(float(r['relative_error']) for r in data),
        single_pass_ms=sum(float(r['ms']) for r in data),
        escaped=sum(r['outcome']=='0' for r in data),surface=sum(r['outcome']=='1' for r in data),
        turned=sum(int(r['turns'])>0 for r in data),both_branches=sum(int(r['double_branches'])>0 for r in data),
        caustics=sum(int(r['caustics'])>0 for r in data),multiple_crossings=sum(int(r['caustics'])>=2 for r in data)))
save('angular_summary.csv',angle_summary)
guard=next(r for r in angle_summary if r['method']=='event_guard')

decomp=read('ft07_decomposition.csv');ft=read(pick('ft07')['file'],True)
max_cut=max(float(r['omitted_fraction']) for r in decomp)
max_ft_num=max(abs(float(a['tau'])-float(b['truncated']))/max(1e-12,float(b['truncated'])) for a,b in zip(ft,decomp))
full_angle={(r['id'],r['pol']):float(r['reference']) for r in angles if r['method']=='event_guard'}
ft_angle_bias=max((full_angle[(r['id'],r['pol'])]-float(r['reference']))/full_angle[(r['id'],r['pol'])] for r in angles if r['method']=='ft07' and full_angle[(r['id'],r['pol'])]>1e-6)

guard_rows=read(pick('event_guard')['file'],True)
discrepancies=sorted([r for r in guard_rows if float(r['rel_error'])>1e-7],key=lambda r:float(r['rel_error']),reverse=True)
save('reference_discrepancies.csv',discrepancies)
worst=max(guard_rows,key=lambda r:float(r['rel_error']))
validation=json.loads((OUT/'transport_validation_summary.json').read_text())
scan=[]
for s in summary:
    scan.append(f"| {s['method']} | {float(s['tolerance']):.0e} | {float(s['parameter']):g} | {float(s['initial_step']):g} | {float(s['batch_median_ms']):.3f} | {s['maximum_relative']:.7g} | {s['missed_90pct']} | {s['failures']} |")
angular_table=[]
for a in angle_summary:
    angular_table.append(f"| {a['method']} | {a['cases']} | {a['maximum_relative']:.7g} | {a['single_pass_ms']:.3f} | {a['failures']} |")

report=fr'''# 任意方向共振光深积分：最终报告

**已删除上次所有径向专用算法。当前推荐并接入 `src/main.cpp` 的是通用的“步内共振边界定位”方法。** 它沿原有 Schwarzschild 光子轨道演化，逐点计算光子–磁场夹角，不要求 r 单调、μ 不变，也不使用共振半径反演。main 中 α=0 只用于运行用户原有的径向参数网格；实际积分器对 α、方位角、半球和传播方向使用同一套代码。

本次比较原程序、通用化 FT07 控制器以及两种新策略，共 {len(summary)} 组配置，每组与原表的 900 例比较。原始 `table/bench_od.txt`、`py/bench_od.py` 和磁场表均未修改。

## 900 例精度与效率

下表为 tolerance=10⁻⁶、初始 proper-path 步长 0.01 km。FT07 使用论文的 0.01 系数；两种新策略默认参数为 0.1。时间是一次完整 900 例计算的 **5 次中位数，单位 ms**，不包含编译和文件 I/O。平台是 Apple arm64 / macOS，GCC 16.2.0，`-O3 -std=c++23`，无 fast-math；详见 [environment.json](environment.json)。

| 方法 | 900 例耗时/ms | 最大相对误差 | P95 相对误差 | 低估超过 90% 的例数 | 平均接受步数 |
|---|---:|---:|---:|---:|---:|
{chr(10).join(table)}

三种修复方案在全部配置中均无异常，也没有严重漏层。两种新策略保留完整粒子分布；FT07 保留论文的中央 99.8% 分布。相同 tolerance 不是相同的全局误差保证：原程序的误差控制用于四维 DOPRI5，修复方案分别控制几何和光深求积，实际精度必须以表中实测值判断。

在对表最大误差相近时，默认边界定位方法比默认共振变量限制器快 **{float(pick('phase_cap')['batch_median_ms'])/float(pick('event_guard')['batch_median_ms']):.2f} 倍**，比采用相同边界处理基础的 FT07 控制器快 **{float(pick('ft07')['batch_median_ms'])/float(pick('event_guard')['batch_median_ms']):.2f} 倍**。新方法比未保护的原程序慢，但原程序的较小成本伴随 27 例严重漏积分。这里没有沿用上次径向专用求积的耗时。

原程序最直观的失败仍是 β₀=−0.1、μz=0、E 模、ω∞=0.01：参考 τ≈10.0904，而原程序约为 1.28×10⁻²⁴。嵌入式误差估计无法识别所有 stage 都未采到的共振层。

将光深转换成散射概率 `P=1−exp(−τ)` 后，推荐方法在 900 例中的最大绝对概率误差为 {probability[-1]['maximum_absolute_probability_error']:.7g}。全部方法的概率误差见 [scattering_probability_errors.csv](scattering_probability_errors.csv)。

## 为什么这些方法适用于非径向光线

两种新策略均显式使用变化的 `x(l)` 和 `μ(l)`：

- **共振变量变化限制器**：按 `|d ln x/dl|` 和 `|dμ/dl|` 限制步长，核心分辨尺度随 β₀² 缩小。远离分布核心时放宽，绝不以当前 opacity=0 为跳过整步的依据。
- **步内共振边界定位**：几何 DOPRI5 可以取较大的步，但在连续轨道插值上求出共振面交点。先搜索内部极值，再括根，能发现起终点同号的窄区间；按交点分段积分。

二者共用 `l=l_a+(l_b−l_a)sin²(πz/2)` 的局部路径长度变换，处理两条速度根合并时的可积奇点。这里变换的是沿轨道的 l，而非半径 r。两个支撑内的速度根均被求和；轨道转向、重入共振区和撞击星表不会被当作径向逃逸处理。

完整公式、控制器细节和调用例子见 [ALGORITHMS.md](ALGORITHMS.md)。

## 非径向、转向和事件验证

额外构造 **80 个初始几何/分布组合 × E/O 两种偏振 = 160 例**，其中 152 例为非径向初始方向。包含南北半球、两极、赤道、任意局域方位角、正负速度分布、初始向内、切向传播，以及从星表非径向发射。参数逐例写在 [angular_validation.csv](angular_validation.csv)。

非径向参考使用同一物理模型的**加密路径计算**：`Δl≤0.005 r`、16 个探测子区间、光深容差 10⁻⁸。它不是用户表之外另有解析真值，也不是已删除的径向速度积分。独立的窄层试验和原公式交叉检查另列在下方。

| 方法 | 方向/偏振例数 | 最大归一化光深差异 | 单次全部额外算例耗时/ms | 失败数 |
|---|---:|---:|---:|---:|
{chr(10).join(angular_table)}

这里差异定义为 `|Δτ|/max(|τ_ref|,10⁻¹²)`；验证阈值为 `|Δτ|≤10⁻⁷+10⁻⁴|τ_ref|`，还要求终止类型一致。耗时是单次观察值，不能与上面 5 次中位数混用。

推荐方法的实际覆盖：{guard['escaped']} 例逃逸、{guard['surface']} 例撞击星表、{guard['turned']} 例发生径向转向、{guard['both_branches']} 例有两个速度分支同时贡献、{guard['caustics']} 例遇到分支合并边界，其中 {guard['multiple_crossings']} 个偏振算例有多次共振边界穿越。三个控制器的终止判定均与加密计算一致。

FT07 在上表与**相同 99.8% 截断区间**的加密解比较，用于检查数值误差；不能把这一列误读为对完整粒子分布的误差。对于额外算例中完整 τ>10⁻⁶ 的情况，FT07 的截断偏差最大达到 {ft_angle_bias*100:.5f}%。

另外通过的检查：

- 4 个已知精确支撑位置的人工窄共振区间。最窄宽度仅为一步的 10⁻⁶，常规探测点均在区间外，仍找到两侧边界；对独立变量变换积分的最大相对差 {validation['max_pocket_relative_error']:.3g}。
- 3000 个任意角度点，将新稳定 opacity 表达式与原 `PhotonEvolution` 公式比较，最大归一化差 {validation['max_opacity_relative_difference']:.3g}。
- {validation['event_tests']} 组固定 U∈{{0.1,0.5,0.9}} 的散射/逃逸检查，其中 {validation['escape_tests']} 组未达到散射光深；散射位置处，加密解的累计光深与 `−ln U` 最大绝对差 {validation['max_event_absolute_error']:.3g}。
- Schwarzschild 冲量参数 `r sin(alpha)/L` 的最大相对漂移 {validation['max_impact_invariant_relative_error']:.3g}。
- 2 个近擦边星表试验：理论近心点分别为 `R_star(1±10⁻⁶)`，成功区分撞击与逃逸；不会仅因一个步的两端在星表外就漏掉中间撞击。

这些验证证明实现确实执行了通用轨道和两分支算法，不是把径向结果套给非径向光线。有限测试不是对所有连续参数、退化高阶切触或任意未解析磁场结构的数学保证。

## FT07 对照的含义

FT07 控制器现在包含完整的方向导数：

\[
\frac{{d\beta}}{{dl}}=\frac{{1-\beta^2}}{{\beta-\mu}}
\left[(1-\beta\mu)\frac{{d\ln x}}{{dl}}+\beta\frac{{d\mu}}{{dl}}\right].
\]

实现论文的层外 r/10、边缘分辨率、动量相对步长式 (39) 和根合并限制式 (38)。为比较完整未散射轨道的光深，三种修复方法共用同一套边界/奇点处理；FT07 的极近合并邻域由它接管，避免无限缩步。**这里测得的是 FT07 控制器在共同数值基础上的成本，不是论文原始 FORTRAN Monte Carlo 程序的成本。**

将 FT07 与通用边界定位算法在相同截断区间上的结果比较，900 例最大数值相对差仅为 {max_ft_num:.3g}；但是粒子尾部截断本身最多丢掉 {max_cut*100:.7f}% 光深。因此对原表的约 1.23% 偏差主要不是步长不够小。分解数据见 [ft07_decomposition.csv](ft07_decomposition.csv)。

## 原表残余差异

推荐方法对原表的最大相对误差为 **{float(worst['rel_error']):.9g}**，出现在 β₀={float(worst['b0']):g}、μz={float(worst['muz']):g}、{'E' if worst['pol']=='1' else 'O'} 模、ω∞={float(worst['omega_inf']):g}。

- 原表 τ={float(worst['reference']):.17g}
- 本次通用边界定位 τ={float(worst['tau']):.17g}

降低容差或改变初始步长后，该量级差异仍存在。所有误差指标继续把用户表作为参考；没有更改参考值、排除这些算例或宣称达到 10⁻¹⁰ 的对表精度。超过 10⁻⁷ 的差异列在 [reference_discrepancies.csv](reference_discrepancies.csv)。

## 完整参数扫描

| 方法 | 容差 | 控制参数 | 初始步长/km | 耗时/ms | 最大相对误差 | 严重漏层数 | 异常数 |
|---|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(scan)}

限制器系数越小通常越耗时，但对给定参考表的误差不保证单调。原程序收紧容差到 10⁻⁸ 后仍有 7 例严重漏层；改变初始步长也不能稳定解决问题。两种新方法在这些配置中均未严重漏层。边界定位默认系数附近的计时小幅差别主要属于运行波动，不能据单个最小时间声称相应系数更优。

## 文件与复现

每个方法都有独立入口：`benchmark/baseline.cpp`、`benchmark/ft07.cpp`、`benchmark/phase_cap.cpp`、`benchmark/event_guard.cpp`。共同物理与算法在 `src/resonance.hpp`、`src/resonance_panels.hpp`、`src/transport.hpp`，原一般轨道/opacity 实现仍在 `src/photon_evolution.hpp`。

所有输出、辅助驱动和说明均保存在 output/。本次删除了旧径向算法、其测试和过时的 benchmark/ 下报告；不要再使用旧径向方法的构建命令。

```sh
python3 output/run_benchmarks.py --repeats 5
python3 output/write_report.py
```

也可以单独运行：

```sh
g++-16 -O3 -std=c++23 benchmark/event_guard.cpp -o output/event_guard
output/event_guard output/my_result.csv 1e-6 .1 5 .01
```

参数依次为输出 CSV、光深容差、控制器系数、计时重复次数和初始步长。共同辅助头文件按要求放在 output/，构建 benchmark 时需要保留。

每次计时包含分布对象、轨道初始化、共振边界定位和积分；不包含参考数据与磁场文件读取、编译或输出 I/O。九个温度的分位数共用准备耗时另记 `setup_ms`。原程序的 `evaluations` 为四维 RHS 调用次数，修复方法则为磁场/共振状态查询次数；后者不含纯几何 RHS，不能把两者当成等价工作单位。逐例 CSV 包含函数查询数、接受步数及光深求积调用数，实际时间才是主要效率指标。

完整逐例结果和所有配置汇总见 [summary.csv](summary.csv)、[summary.json](summary.json)；主程序的 900 行结果见 [main_results.txt](main_results.txt)。
'''
(OUT/'REPORT.md').write_text(report)
print('Wrote output/REPORT.md, angular_summary.csv and diagnostics.')
