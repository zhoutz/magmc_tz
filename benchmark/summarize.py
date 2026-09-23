#!/usr/bin/env python3
"""Summarize transport accuracy and measured per-ray timings (standard library only)."""
import csv
import math
import statistics
import sys
from pathlib import Path

source = Path(sys.argv[1] if len(sys.argv) > 1 else 'benchmark/results.csv')
rows = list(csv.DictReader(source.open()))
methods = list(dict.fromkeys(r['method'] for r in rows))

def percentile(values, p):
    values = sorted(values)
    t = p * (len(values) - 1)
    lo = int(t)
    return values[lo] + (values[min(lo + 1, len(values)-1)]-values[lo])*(t-lo)

lines = [f'数据文件：`{source}`；{len(rows) // len(methods)} 条光线；时间为每例重复运行的中位数。', '']
for label, filtered in [('全部', rows), ('独立径向参考', [r for r in rows if r['radial'] == '1']),
                        ('非径向收敛参考', [r for r in rows if r['radial'] == '0'])]:
    lines += [f'**{label}**', '', '| 方法 | 失败 | 误差>0.1% | 中位相对误差 | 最大相对误差 | 最大概率绝对误差 | 中位耗时 μs | P90耗时 μs | 中位场计算次数 |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for method in methods:
        selected = [r for r in filtered if r['method'] == method]
        good = [r for r in selected if r['status'] == 'ok' and math.isfinite(float(r['relative_error']))]
        if not good:
            continue
        errors = [float(r['relative_error']) for r in good]
        probs = [float(r['probability_error']) for r in good]
        times = [float(r['median_us']) for r in good]
        fields = [int(r['field_calls']) for r in good]
        lines.append(f'| {method} | {len(selected)-len(good)}/{len(selected)} | {sum(e>1e-3 for e in errors)} | {statistics.median(errors):.3g} | {max(errors):.3g} | {max(probs):.3g} | {statistics.median(times):.1f} | {percentile(times,.9):.1f} | {statistics.median(fields):.0f} |')
    lines.append('')
lines += ['**最大误差案例（各方法）**', '', '| 方法 | 案例 | 参考光深 | 算得光深 | 相对误差 |', '|---|---|---:|---:|---:|']
for method in methods:
    good = [r for r in rows if r['method'] == method and r['status']=='ok' and math.isfinite(float(r['relative_error']))]
    if good:
        worst = max(good, key=lambda r: float(r['relative_error']))
        lines.append(f"| {method} | {worst['case']} | {float(worst['tau_reference']):.12g} | {float(worst['tau']):.12g} | {float(worst['relative_error']):.5g} |")
lines += ['', '**历史审查参考案例（本版 RK 在此例未完全漏层）**', '', '| 偏振 | 方法 | 参考光深 | 算得光深 | 耗时 μs |', '|---|---|---:|---:|---:|']
for row in rows:
    if row['case'].startswith('review_'):
        lines.append(f"| {row['pol']} | {row['method']} | {float(row['tau_reference']):.12g} | {float(row['tau']):.12g} | {float(row['median_us']):.1f} |")
failed = [r for r in rows if r['status'] != 'ok']
if failed:
    lines += ['', '**未完成计算（不计入成功计时中位数）**', '']
    lines += [f"- {r['method']} / {r['case']}: {r['status']}" for r in failed]
nonrad = [r for r in rows if r['radial']=='0' and r['method']==methods[0]]
if nonrad:
    lines += ['', f"非径向参考两次加密的最大绝对差：{max(float(r['reference_uncertainty']) for r in nonrad):.6g}。"]
text = '\n'.join(lines) + '\n'
source.with_suffix('.md').write_text(text)
print(text)

# The experimental RK baselines are allowed to fail: exposing those failures
# is the purpose of the comparison. Both event configurations must pass.
event_rows = [r for r in rows if r['method'] in ('event', 'event_tight', 'fast')]
bad = [r for r in event_rows if r['status'] != 'ok' or
       not math.isfinite(float(r['relative_error'])) or float(r['relative_error']) > 1e-6]
if bad:
    raise SystemExit('Event transport regression: ' + ', '.join(r['case']+'/'+r['method'] for r in bad))
