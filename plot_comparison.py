#!/usr/bin/env python3
"""
Generate comparison plots for resonance layer integration methods
"""
import csv
import matplotlib.pyplot as plt
import numpy as np

# Read results
results = []
with open('benchmark_comprehensive.txt', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        results.append({
            'method': row['method'],
            'test_case': row['test_case'],
            'final_tau': float(row['final_tau']),
            'num_steps': int(row['num_steps']),
            'time_ms': float(row['time_ms']),
        })

# Get reference values
reference_values = {}
for r in results:
    if r['method'] == 'Reference':
        reference_values[r['test_case']] = r['final_tau']

# Calculate errors for each method
methods_data = {}
for method in set(r['method'] for r in results):
    if method == 'Reference':
        continue

    method_results = [r for r in results if r['method'] == method]
    errors = []
    steps = []
    times = []

    for r in method_results:
        ref = reference_values[r['test_case']]
        err = abs(r['final_tau'] - ref) / abs(ref) * 100
        errors.append(err)
        steps.append(r['num_steps'])
        times.append(r['time_ms'])

    methods_data[method] = {
        'mean_error': np.mean(errors),
        'max_error': np.max(errors),
        'mean_steps': np.mean(steps),
        'mean_time': np.mean(times)
    }

# Create figure with subplots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Resonance Layer Integration Methods Comparison', fontsize=16, fontweight='bold')

# Plot 1: Error vs Steps (log scale)
ax1 = axes[0, 0]
for method, data in methods_data.items():
    color = 'red' if 'Original' in method else 'blue' if 'DtauLimit' in method else 'gray'
    marker = 'o' if 'DtauLimit' in method else 's' if 'Original' in method else '^'
    size = 100 if 'DtauLimit_0.01' in method else 60
    ax1.scatter(data['mean_steps'], data['mean_error'],
               s=size, alpha=0.7, color=color, marker=marker, label=method)

ax1.set_xlabel('Mean Steps', fontsize=12)
ax1.set_ylabel('Mean Relative Error (%)', fontsize=12)
ax1.set_title('Accuracy vs Computational Cost')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=8, loc='best')

# Plot 2: Error comparison (bar chart)
ax2 = axes[0, 1]
method_names = list(methods_data.keys())
mean_errors = [methods_data[m]['mean_error'] for m in method_names]
colors = ['red' if 'Original' in m else 'green' if 'DtauLimit_0.01' in m else 'blue' if 'DtauLimit' in m else 'gray'
          for m in method_names]

bars = ax2.barh(range(len(method_names)), mean_errors, color=colors, alpha=0.7)
ax2.set_yticks(range(len(method_names)))
ax2.set_yticklabels(method_names, fontsize=9)
ax2.set_xlabel('Mean Relative Error (%)', fontsize=12)
ax2.set_title('Accuracy Comparison')
ax2.grid(True, alpha=0.3, axis='x')
ax2.axvline(x=0.01, color='green', linestyle='--', alpha=0.5, label='Target: 0.01%')
ax2.legend()

# Plot 3: Time comparison (bar chart)
ax3 = axes[1, 0]
mean_times = [methods_data[m]['mean_time'] for m in method_names]
bars = ax3.barh(range(len(method_names)), mean_times, color=colors, alpha=0.7)
ax3.set_yticks(range(len(method_names)))
ax3.set_yticklabels(method_names, fontsize=9)
ax3.set_xlabel('Mean Time (ms)', fontsize=12)
ax3.set_title('Computational Time Comparison')
ax3.set_xscale('log')
ax3.grid(True, alpha=0.3, axis='x')

# Plot 4: Efficiency metric (error × time)
ax4 = axes[1, 1]
efficiency = [methods_data[m]['mean_error'] * methods_data[m]['mean_time'] for m in method_names]
bars = ax4.barh(range(len(method_names)), efficiency, color=colors, alpha=0.7)
ax4.set_yticks(range(len(method_names)))
ax4.set_yticklabels(method_names, fontsize=9)
ax4.set_xlabel('Efficiency Metric (Error × Time)', fontsize=12)
ax4.set_title('Overall Efficiency (Lower is Better)')
ax4.set_xscale('log')
ax4.grid(True, alpha=0.3, axis='x')

# Highlight best method
best_method = min(method_names, key=lambda m: methods_data[m]['mean_error'] * methods_data[m]['mean_time'])
best_idx = method_names.index(best_method)
for ax in [ax2, ax3, ax4]:
    ax.get_children()[best_idx].set_edgecolor('green')
    ax.get_children()[best_idx].set_linewidth(3)

plt.tight_layout()
plt.savefig('resonance_methods_comparison.png', dpi=300, bbox_inches='tight')
print("Plot saved to: resonance_methods_comparison.png")

# Print summary
print("\n" + "="*80)
print("BEST METHOD: " + best_method)
print("="*80)
print(f"Mean Error: {methods_data[best_method]['mean_error']:.6f}%")
print(f"Mean Steps: {methods_data[best_method]['mean_steps']:.0f}")
print(f"Mean Time:  {methods_data[best_method]['mean_time']:.2f} ms")
print(f"Efficiency: {methods_data[best_method]['mean_error'] * methods_data[best_method]['mean_time']:.6f}")
