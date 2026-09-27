#!/usr/bin/env python3
"""
Analysis script for polarization output from fd11.cpp
Plots the Stokes parameters and polarization degree
"""

import numpy as np
import matplotlib.pyplot as plt

# Read data
data = np.loadtxt('output/fd11.txt')
omega_inf = data[:, 0]
muk = data[:, 1]
Q = data[:, 2]
U = data[:, 3]
V = data[:, 4]

# Compute polarization degree and angle
P_linear = np.sqrt(Q**2 + U**2)
P_total = np.sqrt(Q**2 + U**2 + V**2)
chi = 0.5 * np.arctan2(U, Q) * 180 / np.pi  # Polarization angle in degrees

# Create plots
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Plot 1: Q vs U (Poincaré disk projection)
ax = axes[0, 0]
ax.scatter(Q, U, c=omega_inf, s=5, alpha=0.5, cmap='viridis')
ax.set_xlabel('Q')
ax.set_ylabel('U')
ax.set_title('Q-U Plane')
ax.axis('equal')
ax.grid(True, alpha=0.3)
circle = plt.Circle((0, 0), 1, fill=False, color='red', linestyle='--', alpha=0.5)
ax.add_patch(circle)

# Plot 2: V distribution
ax = axes[0, 1]
ax.hist(V, bins=50, alpha=0.7, edgecolor='black')
ax.set_xlabel('V (Circular Polarization)')
ax.set_ylabel('Count')
ax.set_title('V Distribution')
ax.grid(True, alpha=0.3)

# Plot 3: Linear polarization degree
ax = axes[0, 2]
ax.hist(P_linear, bins=50, alpha=0.7, edgecolor='black', color='green')
ax.set_xlabel('Linear Polarization Degree')
ax.set_ylabel('Count')
ax.set_title('P_linear = sqrt(Q² + U²)')
ax.grid(True, alpha=0.3)

# Plot 4: Total polarization degree
ax = axes[1, 0]
ax.hist(P_total, bins=50, alpha=0.7, edgecolor='black', color='orange')
ax.set_xlabel('Total Polarization Degree')
ax.set_ylabel('Count')
ax.set_title('P_total = sqrt(Q² + U² + V²)')
ax.grid(True, alpha=0.3)

# Plot 5: Polarization vs energy
ax = axes[1, 1]
ax.scatter(omega_inf, P_total, c=np.abs(muk), s=5, alpha=0.5, cmap='plasma')
ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Total Polarization')
ax.set_title('Polarization vs Energy')
ax.set_xscale('log')
ax.grid(True, alpha=0.3)
cbar = plt.colorbar(ax.collections[0], ax=ax)
cbar.set_label('|cos θ|')

# Plot 6: Polarization vs viewing angle
ax = axes[1, 2]
ax.scatter(muk, P_total, c=omega_inf, s=5, alpha=0.5, cmap='viridis')
ax.set_xlabel('cos θ_view')
ax.set_ylabel('Total Polarization')
ax.set_title('Polarization vs Viewing Angle')
ax.grid(True, alpha=0.3)
cbar = plt.colorbar(ax.collections[0], ax=ax)
cbar.set_label('Energy (keV)')

plt.tight_layout()
plt.savefig('output/fd11_polarization_analysis.png', dpi=150)
print('Saved plot to output/fd11_polarization_analysis.png')

# Print statistics
print('\n=== Polarization Statistics ===')
print(f'Total photons: {len(Q)}')
print(f'\nMean Q: {np.mean(Q):.4f} ± {np.std(Q):.4f}')
print(f'Mean U: {np.mean(U):.4f} ± {np.std(U):.4f}')
print(f'Mean V: {np.mean(V):.4f} ± {np.std(V):.4f}')
print(f'\nMean linear polarization: {np.mean(P_linear):.4f} ± {np.std(P_linear):.4f}')
print(f'Mean total polarization: {np.mean(P_total):.4f} ± {np.std(P_total):.4f}')
print(f'\nMedian linear polarization: {np.median(P_linear):.4f}')
print(f'Median total polarization: {np.median(P_total):.4f}')

# Check for photons with high circular polarization
high_V = np.abs(V) > 0.1
if np.any(high_V):
    print(f'\nPhotons with |V| > 0.1: {np.sum(high_V)} ({100*np.sum(high_V)/len(V):.2f}%)')
    print(f'Max |V|: {np.max(np.abs(V)):.4f}')
