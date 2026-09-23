#!/usr/bin/env python3
"""Paired speed and accuracy summary for the previous and fast transports."""
import csv
import math
import statistics
from pathlib import Path

root = Path(__file__).resolve().parent
rows = list(csv.DictReader((root/'fast_results.csv').open()))
old = {r['case']: r for r in rows if r['method']=='event'}
fast = [r for r in rows if r['method']=='fast']
assert len(old)==len(fast) and fast
bad = [r for r in fast if r['status']!='ok' or not math.isfinite(float(r['relative_error']))
       or float(r['relative_error'])>1e-6]
if bad:
    raise SystemExit('Fast regression: '+', '.join(r['case'] for r in bad))
lines=['# 加速实测汇总','',f'{len(fast)} 条相同初始光线；新旧方案均单线程。原始数据：[fast_results.csv](fast_results.csv)。','',
       '| 样本 | 旧方案中位耗时 μs | 新方案中位耗时 μs | 逐例加速比中位数 | 总耗时比 | 新方案最大光深相对误差 |',
       '|---|---:|---:|---:|---:|---:|']
for label, selected in [('全部',fast),('径向',[r for r in fast if r['radial']=='1']),
                        ('非径向',[r for r in fast if r['radial']=='0'])]:
    if not selected: continue
    olds=[float(old[r['case']]['median_us']) for r in selected]
    news=[float(r['median_us']) for r in selected]
    ratios=[a/b for a,b in zip(olds,news)]
    lines.append(f'| {label} | {statistics.median(olds):.1f} | {statistics.median(news):.1f} | '
                 f'{statistics.median(ratios):.2f}× | {sum(olds)/sum(news):.2f}× | '
                 f'{max(float(r["relative_error"]) for r in selected):.3g} |')
lines+=['','总耗时比为同组所有案例耗时中位数之和的比值；不把“中位耗时之比”误称为逐例加速比。',
        '径向使用独立速度积分参考；非径向使用旧方案加密后的收敛参考。',
        f'新方案最大散射概率绝对误差：{max(float(r["probability_error"]) for r in fast):.3g}。']
ensemble = root/'ensemble_speed.csv'
if ensemble.exists():
    data=list(csv.DictReader(ensemble.open()))
    lines+=['','## 完整多次散射程序','',
            '| beta0 | 初始偏振 | 光子数 | 旧方案中位耗时 s | 新方案中位耗时 s | 加速 | 终止/散射计数一致 |',
            '|---:|---|---:|---:|---:|---:|---|']
    sumold=sumfast=0
    for key in dict.fromkeys((r['beta0'],r['mode']) for r in data):
        selected=[r for r in data if (r['beta0'],r['mode'])==key]
        a=[r for r in selected if r['method']=='event'];b=[r for r in selected if r['method']=='fast']
        ta=statistics.median(float(r['seconds']) for r in a)
        tb=statistics.median(float(r['seconds']) for r in b)
        fields=['escaped','absorbed','scatterings']
        same=all(tuple(r[k] for k in fields)==tuple(a[0][k] for k in fields) for r in selected)
        lines.append(f'| {key[0]} | {key[1]} | {a[0]["photons"]} | {ta:.6f} | {tb:.6f} | {ta/tb:.2f}× | {"是" if same else "否"} |')
        sumold+=ta;sumfast+=tb
    lines += ['', f'四组总耗时比：{sumold/sumfast:.2f}×。包含进程启动、磁场表读取和输出；各组新旧方法交替运行。',
              '计数一致只是附加检查，精度由独立光深参考与逐路径/散射位置回归判断。']
text='\n'.join(lines)+'\n'
(root/'fast_results.md').write_text(text)
print(text)
