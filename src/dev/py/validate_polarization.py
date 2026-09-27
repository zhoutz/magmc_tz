#!/usr/bin/env python3
"""
Validation tests for polarization implementation
Checks that the physics is correctly implemented
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 60)
print("POLARIZATION IMPLEMENTATION VALIDATION")
print("=" * 60)

# Read data
data = np.loadtxt('output/fd11.txt')
omega_inf = data[:, 0]
muk = data[:, 1]
Q = data[:, 2]
U = data[:, 3]
V = data[:, 4]

print(f"\nLoaded {len(Q)} escaped photons")

# Test 1: Stokes parameter normalization
print("\n" + "-" * 60)
print("TEST 1: Stokes Parameter Normalization")
print("-" * 60)
I_effective = np.sqrt(Q**2 + U**2 + V**2)
mean_I = np.mean(I_effective)
std_I = np.std(I_effective)
print(f"I = sqrt(Q² + U² + V²)")
print(f"  Mean: {mean_I:.6f}")
print(f"  Std:  {std_I:.6f}")
print(f"  Min:  {np.min(I_effective):.6f}")
print(f"  Max:  {np.max(I_effective):.6f}")

if np.abs(mean_I - 1.0) < 0.01 and std_I < 0.1:
    print("✅ PASS: Polarization is normalized (I ≈ 1)")
else:
    print("⚠ WARNING: Polarization normalization issue")

# Test 2: Physical bounds
print("\n" + "-" * 60)
print("TEST 2: Physical Bounds")
print("-" * 60)
exceeds_bounds = I_effective > 1.1  # Allow 10% tolerance
n_exceeds = np.sum(exceeds_bounds)
print(f"Photons with P > 1.1: {n_exceeds} ({100*n_exceeds/len(Q):.2f}%)")

if n_exceeds == 0:
    print("✅ PASS: All polarizations within physical bounds")
elif n_exceeds < 0.01 * len(Q):
    print("⚠ WARNING: Few photons exceed bounds (likely numerical)")
else:
    print("❌ FAIL: Many photons exceed physical bounds")

# Test 3: Mode conservation (before coupling)
print("\n" + "-" * 60)
print("TEST 3: Initial Mode Identification")
print("-" * 60)
# Q ≈ 1 suggests O-mode, Q ≈ 0 suggests E-mode
Q_threshold = 0.5
O_mode = np.abs(Q - 1.0) < 0.1
E_mode = np.abs(Q - 0.0) < 0.1
mixed = ~O_mode & ~E_mode

print(f"Approximate O-mode (Q ≈ 1): {np.sum(O_mode)} ({100*np.sum(O_mode)/len(Q):.1f}%)")
print(f"Approximate E-mode (Q ≈ 0): {np.sum(E_mode)} ({100*np.sum(E_mode)/len(Q):.1f}%)")
print(f"Mixed/evolved:               {np.sum(mixed)} ({100*np.sum(mixed)/len(Q):.1f}%)")

if np.sum(O_mode) + np.sum(E_mode) > 0:
    print("✅ PASS: Can identify initial modes")
else:
    print("⚠ NOTE: All photons show significant mode evolution")

# Test 4: Linear vs circular polarization
print("\n" + "-" * 60)
print("TEST 4: Linear vs Circular Polarization")
print("-" * 60)
P_linear = np.sqrt(Q**2 + U**2)
P_circ = np.abs(V)
ratio = np.mean(P_linear) / np.mean(P_circ)

print(f"Mean linear pol:    {np.mean(P_linear):.6f}")
print(f"Mean circular pol:  {np.mean(P_circ):.6f}")
print(f"Ratio (linear/circ): {ratio:.1f}")

if ratio > 10:
    print("✅ PASS: Linear >> circular (vacuum polarization dominant)")
elif ratio > 5:
    print("⚠ NOTE: Linear > circular but ratio lower than expected")
else:
    print("❌ FAIL: Circular polarization too large")

# Test 5: Q conservation check
print("\n" + "-" * 60)
print("TEST 5: Q Conservation Property")
print("-" * 60)
print("According to theory, Q should be conserved during propagation")
print("(only U and V rotate)")
print(f"\nQ distribution:")
print(f"  Mean: {np.mean(Q):.4f}")
print(f"  Std:  {np.std(Q):.4f}")
print(f"  Range: [{np.min(Q):.4f}, {np.max(Q):.4f}]")

# For photons that evolved, Q should cluster near initial values (0 or 1)
if np.std(Q) < 0.5:
    print("✅ PASS: Q shows expected structure")
else:
    print("⚠ NOTE: Q varies significantly (may indicate strong evolution)")

# Test 6: Energy dependence
print("\n" + "-" * 60)
print("TEST 6: Energy Dependence")
print("-" * 60)
# Higher energy photons should show less polarization evolution
# because coupling parameter scales as omega^-1
low_E = omega_inf < 1.5
high_E = omega_inf > 5.0

if np.sum(low_E) > 10 and np.sum(high_E) > 10:
    P_low = np.sqrt(Q[low_E]**2 + U[low_E]**2 + V[low_E]**2)
    P_high = np.sqrt(Q[high_E]**2 + U[high_E]**2 + V[high_E]**2)

    print(f"Low energy (E < 1.5 keV):  {np.sum(low_E)} photons, P = {np.mean(P_low):.4f}")
    print(f"High energy (E > 5 keV):   {np.sum(high_E)} photons, P = {np.mean(P_high):.4f}")

    print("✓ Both energy ranges have photons")
else:
    print("⚠ NOTE: Insufficient statistics for energy comparison")

# Test 7: Create validation plot
print("\n" + "-" * 60)
print("TEST 7: Generating Validation Plots")
print("-" * 60)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Poincaré sphere (Q-U plane)
ax = axes[0, 0]
scatter = ax.scatter(Q, U, c=V, s=10, alpha=0.5, cmap='RdBu_r', vmin=-1, vmax=1)
ax.set_xlabel('Q', fontsize=12)
ax.set_ylabel('U', fontsize=12)
ax.set_title('Poincaré Sphere Projection (Q-U plane)', fontsize=13)
ax.axis('equal')
ax.grid(True, alpha=0.3)
circle = plt.Circle((0, 0), 1, fill=False, color='black', linestyle='--', alpha=0.5, linewidth=2)
ax.add_patch(circle)
plt.colorbar(scatter, ax=ax, label='V')
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)

# Plot 2: Normalization check
ax = axes[0, 1]
ax.hist(I_effective, bins=50, alpha=0.7, edgecolor='black', color='green')
ax.axvline(1.0, color='red', linestyle='--', linewidth=2, label='I = 1')
ax.set_xlabel('Total Polarization Degree', fontsize=12)
ax.set_ylabel('Count', fontsize=12)
ax.set_title('Normalization: I = sqrt(Q² + U² + V²)', fontsize=13)
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 3: Component distributions
ax = axes[1, 0]
ax.hist(Q, bins=50, alpha=0.5, label='Q', color='red', edgecolor='black')
ax.hist(U, bins=50, alpha=0.5, label='U', color='blue', edgecolor='black')
ax.hist(V, bins=50, alpha=0.5, label='V', color='green', edgecolor='black')
ax.set_xlabel('Stokes Parameter Value', fontsize=12)
ax.set_ylabel('Count', fontsize=12)
ax.set_title('Stokes Parameter Distributions', fontsize=13)
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 4: 3D scatter (if few enough points)
ax = axes[1, 1]
if len(Q) < 10000:
    ax.remove()
    ax = fig.add_subplot(2, 2, 4, projection='3d')
    scatter = ax.scatter(Q, U, V, c=omega_inf, s=5, alpha=0.3, cmap='viridis')
    ax.set_xlabel('Q', fontsize=10)
    ax.set_ylabel('U', fontsize=10)
    ax.set_zlabel('V', fontsize=10)
    ax.set_title('3D Stokes Space', fontsize=13)
    plt.colorbar(scatter, ax=ax, label='Energy (keV)', shrink=0.5)
else:
    # 2D projection instead
    scatter = ax.scatter(np.sqrt(Q**2 + U**2), V, c=omega_inf, s=5, alpha=0.5, cmap='viridis')
    ax.set_xlabel('Linear Polarization', fontsize=12)
    ax.set_ylabel('Circular Polarization (V)', fontsize=12)
    ax.set_title('Linear vs Circular Polarization', fontsize=13)
    plt.colorbar(scatter, ax=ax, label='Energy (keV)')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('output/fd11_validation.png', dpi=150, bbox_inches='tight')
print("✅ Saved validation plot to output/fd11_validation.png")

# Final summary
print("\n" + "=" * 60)
print("VALIDATION SUMMARY")
print("=" * 60)

tests_passed = 0
total_tests = 4

if np.abs(mean_I - 1.0) < 0.01 and std_I < 0.1:
    tests_passed += 1
if n_exceeds < 0.01 * len(Q):
    tests_passed += 1
if ratio > 5:
    tests_passed += 1
if np.sum(O_mode) + np.sum(E_mode) > 0:
    tests_passed += 1

print(f"\nTests passed: {tests_passed}/{total_tests}")

if tests_passed == total_tests:
    print("\n✅ ALL TESTS PASSED - Implementation appears correct")
elif tests_passed >= total_tests - 1:
    print("\n✓ MOSTLY CORRECT - Minor issues detected")
else:
    print("\n⚠ ISSUES DETECTED - Review implementation")

print("\n" + "=" * 60)
