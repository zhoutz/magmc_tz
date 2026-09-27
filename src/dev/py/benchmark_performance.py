#!/usr/bin/env python3
"""
Performance comparison between original main.cpp and fd11.cpp with polarization
"""

import subprocess
import time
import os

def compile_code(source, output, desc):
    """Compile a C++ source file"""
    cmd = [
        'g++', '-std=c++23', '-O3', '-o', output, source,
        '-I', 'src', '-I/opt/homebrew/include',
        '-L/opt/homebrew/lib', '-lgsl', '-lgslcblas'
    ]
    print(f"Compiling {desc}...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Compilation failed: {result.stderr}")
        return False
    print(f"✓ {desc} compiled successfully")
    return True

def run_benchmark(executable, output_file, n_photons=1000):
    """Run a benchmark and measure execution time"""
    # Modify the code temporarily to use fewer photons
    start = time.time()
    result = subprocess.run([executable], capture_output=True, text=True)
    end = time.time()

    if result.returncode != 0:
        print(f"Execution failed: {result.stderr}")
        return None

    elapsed = end - start

    # Count escaped photons
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            lines = [l for l in f if not l.startswith('#')]
            n_escaped = len(lines)
    else:
        n_escaped = 0

    return elapsed, n_escaped

def main():
    print("=" * 60)
    print("Performance Comparison: main.cpp vs fd11.cpp")
    print("=" * 60)

    # Check if original main binary exists, if not compile it
    if not os.path.exists('build/main'):
        print("\nCompiling original code for comparison...")
        compile_code('src/dev/main.cpp', 'build/main', 'main.cpp')

    print("\nBoth implementations use N=10000 photons (hardcoded)")
    print("\nNote: This is a rough comparison. For precise benchmarking,")
    print("multiple runs with statistics would be needed.")
    print("\nPress Enter to continue...")
    input()

    # Run original code
    print("\n" + "-" * 60)
    print("Running ORIGINAL code (without polarization)...")
    print("-" * 60)
    start_orig = time.time()
    result_orig = subprocess.run(['./build/main'], capture_output=True, text=True)
    end_orig = time.time()
    time_orig = end_orig - start_orig

    if os.path.exists('output/main.txt'):
        with open('output/main.txt', 'r') as f:
            n_orig = len([l for l in f if not l.startswith('#')])
    else:
        n_orig = 0

    print(f"✓ Original code completed in {time_orig:.2f} seconds")
    print(f"  Escaped photons: {n_orig}")

    # Run fd11 code
    print("\n" + "-" * 60)
    print("Running FD11 code (WITH polarization)...")
    print("-" * 60)
    start_fd11 = time.time()
    result_fd11 = subprocess.run(['./build/fd11'], capture_output=True, text=True)
    end_fd11 = time.time()
    time_fd11 = end_fd11 - start_fd11

    if os.path.exists('output/fd11.txt'):
        with open('output/fd11.txt', 'r') as f:
            n_fd11 = len([l for l in f if not l.startswith('#')])
    else:
        n_fd11 = 0

    print(f"✓ FD11 code completed in {time_fd11:.2f} seconds")
    print(f"  Escaped photons: {n_fd11}")

    # Summary
    print("\n" + "=" * 60)
    print("PERFORMANCE SUMMARY")
    print("=" * 60)
    print(f"Original (no polarization):  {time_orig:8.2f} seconds")
    print(f"FD11 (with polarization):    {time_fd11:8.2f} seconds")
    print(f"Overhead:                    {time_fd11 - time_orig:8.2f} seconds ({100*(time_fd11/time_orig - 1):.1f}%)")
    print(f"\nSpeedup factor:              {time_fd11/time_orig:.2f}x")

    if time_fd11 / time_orig < 1.5:
        print("\n✅ EXCELLENT: Overhead is less than 50%")
    elif time_fd11 / time_orig < 2.0:
        print("\n✓ GOOD: Overhead is less than 100%")
    else:
        print("\n⚠ WARNING: Significant overhead detected")

    print("\nNote: The overhead includes:")
    print("  - Computing η_couple at each geodesic step")
    print("  - Polarization ODE integration (when active)")
    print("  - Additional (Q,U,V) storage and output")
    print("=" * 60)

if __name__ == '__main__':
    main()
