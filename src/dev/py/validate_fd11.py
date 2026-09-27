#!/usr/bin/env python3
"""Reproducible FD11 regression/benchmark; run from the repository root.

Source artifacts stay in src/dev; all run data/logs are written to output.
Build both binaries first (see summary/fd11.md).
"""
import json
import math
import pathlib
import re
import statistics
import subprocess
import time

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)


def read_rows(path):
    with path.open() as stream:
        return [list(map(float, line.split())) for line in stream
                if line.strip() and not line.startswith("#")]


def run(label, args):
    path = OUT / (label + ".txt")
    command = [str(ROOT / "build/fd11"), "--output", str(path), *args]
    start = time.perf_counter()
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                               timeout=180)
    elapsed = time.perf_counter() - start
    (OUT / (label + ".log")).write_text(completed.stderr)
    completed.check_returncode()
    rows = read_rows(path)
    assert rows and all(len(row) == 5 for row in rows)
    assert all(all(math.isfinite(v) for v in row) for row in rows)
    assert all(row[0] > 0 and -1 <= row[1] <= 1 for row in rows)
    norm_error = max(abs(sum(v * v for v in row[2:]) - 1) for row in rows)
    assert norm_error < 1e-9, norm_error
    counters = dict(re.findall(r"(\w+)=([\d.]+)", completed.stderr))
    assert int(counters["escaped"]) + int(counters["absorbed"]) == int(counters["photons"])
    assert len(rows) == int(counters["escaped"])
    result = dict(label=label, seconds=elapsed, norm_error=norm_error, counters=counters)
    print(json.dumps(result), flush=True)
    return rows, result


results = []
baseline_times = []
enabled_times = []
disabled_times = []
for repeat in range(3):
    start = time.perf_counter()
    with (OUT / "fd11_original_stdout.log").open("w") as stream:
        subprocess.run([str(ROOT / "build/main_fd11_baseline")], cwd=ROOT,
                       stdout=stream, stderr=subprocess.PIPE, timeout=180, check=True)
    baseline_times.append(time.perf_counter() - start)
    original = read_rows(OUT / "main.txt")
    disabled, record = run(f"fd11_disabled_{repeat}", ["--no-polarization"])
    assert [row[:2] for row in disabled] == original, "Original transport changed"
    disabled_times.append(record["seconds"])
    results.append(record)
    _, record = run(f"fd11_enabled_{repeat}", [])
    enabled_times.append(record["seconds"])
    results.append(record)

for mode in ["E", "O"]:
    for energy in [0.1, 1, 10]:
        _, record = run(f"fd11_sweep_{mode}_{energy}",
                        ["--photons", "1000", "--mode", mode, "--energy", str(energy),
                         "--seed", "4193"])
        results.append(record)

_, record = run("fd11_positrons", ["--photons", "1000", "--charge", "positron", "--seed", "8171"])
results.append(record)

# Same random stream, tolerance refinement. Photon paths may differ when a
# coupled photon scatters; require aggregate finite output rather than silently
# equating rows after such a branching event. Single-ray convergence is tested
# independently by fd11_test.cpp.
_, record = run("fd11_tight", ["--photons", "1000", "--pol-tol", "1e-8"])
results.append(record)

report = dict(original_seconds=baseline_times, enabled_seconds=enabled_times,
              disabled_seconds=disabled_times,
              enabled_over_original=statistics.median(enabled_times) / statistics.median(baseline_times),
              runs=results)
(OUT / "fd11_validation.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps({k: v for k, v in report.items() if k != "runs"}, indent=2))
