"""Run from the repository root with a Python that has NumPy installed.

  /opt/homebrew/Caskroom/miniforge/base/bin/python src/bench/run_bench.py
All timings are sequential, warmed up, and measured inside C++ (table loading
and process startup excluded). No timing results from parallel runs are used.
"""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--repeat", type=int, default=5)
parser.add_argument("--skip-build", action="store_true")
args = parser.parse_args()
methods = ["limited", "event_qags", "event_sqrt", "global_sqrt", "ref_general", "ref_regular"]
env = os.environ.copy()
env.pop("OD_SCALE", None)
env.pop("OD_REF_SCALE", None)
Path("build").mkdir(exist_ok=True)
Path("output").mkdir(exist_ok=True)

def run(command, **kwargs):
    p = subprocess.run(command, text=True, capture_output=True, env=kwargs.pop("env", env), **kwargs)
    if p.returncode:
        raise RuntimeError(f"{command}\n{p.stdout}\n{p.stderr}")
    return p.stdout + p.stderr

if not args.skip_build:
    for name in methods + ["base", "ref", "base_timed", "ref_timed", "validate_od"]:
        source = "original_timing" if name.endswith("_timed") else name
        cmd = ["g++-16", f"src/bench/{source}.cpp", "-o", f"build/{name}", "-std=c++23", "-O3",
               "-lgsl", "-I/opt/homebrew/include", "-L/opt/homebrew/lib"]
        if name == "ref_timed":
            cmd.append("-DBENCH_REFERENCE")
        run(cmd)

# Generate the unmodified originals too, so a clean output/ directory is enough.
run(["build/base"])
run(["build/ref"])

timings = {}
with open("output/timing_runs.txt", "w") as log:
    for suite in ["radial", "nonradial", "stress"]:
        for name in (["base_timed", "ref_timed"] if suite == "radial" else []) + methods:
            command = [f"build/{name}", suite]
            run(command)  # warm-up
            times = []
            for repeat in range(args.repeat):
                result = run(command)
                log.write(result)
                times.append(float(re.search(r"seconds=([0-9.]+)", result).group(1)))
            suffix = "" if suite == "radial" else "_" + suite
            Path(f"output/{name}{suffix}_log.txt").write_text(result)
            timings[f"{name}/{suite}"] = dict(median=statistics.median(times), min=min(times), max=max(times), runs=times)
            print(name, suite, timings[f"{name}/{suite}"]["median"], flush=True)

comparisons = []
for name in ["base", "ref", "base_timed", "ref_timed"] + methods:
    result = run([sys.executable, "py/cmp_od.py", name])
    comparisons.append(result)
Path("output/cmp_od_results.txt").write_text("".join(comparisons))

# Each tight run halves/quarters both the spatial cap and numerical tolerances.
for name, var in [("global_sqrt", "OD_SCALE"), ("ref_regular", "OD_REF_SCALE")]:
    for suite in ["radial", "nonradial", "stress"]:
        result = run([f"build/{name}", suite, "_tight"], env={**env, var: "0.25"})
        Path(f"output/{name}_{suite}_tight_log.txt").write_text(result)

metrics = {}
for original, timed in [("base", "base_timed"), ("ref", "ref_timed")]:
    if not np.array_equal(np.loadtxt(f"output/{original}.txt"), np.loadtxt(f"output/{timed}.txt")):
        raise RuntimeError("timing wrapper changed original results: " + original)
for suite in ["radial", "nonradial", "stress"]:
    suffix = "" if suite == "radial" else "_" + suite
    truth_file = "table/bench_od.txt" if suite == "radial" else f"output/ref_regular{suffix}_tight.txt"
    truth_data = np.loadtxt(truth_file)
    truth = truth_data[:, -1]
    for name in (["base", "ref", "base_timed", "ref_timed"] if suite == "radial" else []) + methods + ["global_sqrt_tight", "ref_regular_tight"]:
        filename = (name[:-6] + suffix + "_tight") if name.endswith("_tight") else name + suffix
        data = np.loadtxt(f"output/{filename}.txt")
        if data.shape != truth_data.shape or not np.allclose(data[:, :-1], truth_data[:, :-1], rtol=0, atol=1e-14):
            raise RuntimeError("case ordering mismatch: " + filename)
        values = data[:, -1]
        if not np.all(np.isfinite(values)):
            raise RuntimeError("nonfinite benchmark: " + filename)
        error = np.abs(values - truth)
        mask = np.abs(truth) > 1e-10
        rel = error[mask] / np.abs(truth[mask])
        metrics[f"{name}/{suite}"] = dict(n=len(values), max_abs=float(max(error)),
            max_rel=float(max(rel)), mean_rel=float(np.mean(rel)), zero_missed=int(np.sum((values == 0) & mask)),
            max_abs_zero_truth=float(np.max(error[~mask], initial=0)),
            worst_rel_case=int(np.flatnonzero(mask)[rel.argmax()]))

Path("output/metrics.txt").write_text(json.dumps(dict(timings=timings, errors=metrics), indent=2) + "\n")
Path("output/validation_log.txt").write_text(run(["build/validate_od"]))
print("Saved output/timing_runs.txt, cmp_od_results.txt, metrics.txt and validation files.")
