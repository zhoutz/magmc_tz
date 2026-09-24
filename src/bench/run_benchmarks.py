"""Run from the repository root; all benchmark code stays in src/bench.

Example: /opt/homebrew/Caskroom/miniforge/base/bin/python src/bench/run_benchmarks.py
Use --analyze to regenerate the Chinese report from saved output only.
"""
from pathlib import Path
import argparse
import subprocess
import sys

METHODS = {"ft07": ".01", "geometric": ".001", "velocity_panels": ".2",
           "base_general": ".2", "ref_general": ".01"}


def command(args, log=None):
    print(" ".join(map(str, args)), flush=True)
    p = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(p.stdout, end="", flush=True)
    if log is not None:
        log.write("$ " + " ".join(map(str, args)) + "\n" + p.stdout)
        log.flush()
    if p.returncode:
        raise RuntimeError(f"Command failed ({p.returncode}): {args}")


def run():
    Path("build").mkdir(exist_ok=True)
    Path("output").mkdir(exist_ok=True)
    with open("output/benchmark_run.txt", "w") as log:
        for name in [*METHODS, "radial_truth", "check_transport", "time_original"]:
            flags = ["g++-16", f"src/bench/{name}.cpp", "-o", f"build/{name}",
                     "-std=c++23", "-O3", "-Wall", "-Wextra", "-Wno-psabi",
                     "-I/opt/homebrew/include", "-L/opt/homebrew/lib", "-lgsl", "-lgslcblas"]
            command(flags, log)
            if name == "time_original":
                flags[flags.index("-o")+1] = "build/time_ref"
                command(flags + ["-DTIME_REF"], log)
        command(["build/time_original"], log)
        command(["build/time_ref"], log)
        command(["build/radial_truth"], log)
        command([sys.executable, "src/bench/audit_table.py"], log)
        for suite in ["radial", "angular", "general"]:
            for name, resolution in METHODS.items():
                stem = name if suite == "radial" else name + "_" + suite
                command(["build/"+name, suite, stem, "1e-6", resolution, "3"], log)
            for name, stem, tol, step, repeats, probes in [
                    ("ref_general", "ref_tight", "1e-8", ".005", "3", "16"),
                    ("velocity_panels", "panels_tight", "1e-8", ".1", "1", "16"),
                    ("ft07", "ft07_refined", "1e-8", ".005", "1", "16")]:
                command(["build/"+name, suite, stem+"_"+suite, tol, step, repeats, ".01", probes], log)
        for suite in ["angular", "general"]:
            for name, stem, tol, step in [("ref_general", "ref_finer", "3e-9", ".0025"),
                                          ("velocity_panels", "panels_finer", "1e-9", ".05")]:
                command(["build/"+name, suite, stem+"_"+suite, tol, step, "1", ".01", "32"], log)
        # Initial-step and event-probe sensitivity, with all other settings fixed.
        command(["build/velocity_panels", "angular", "panels_initial", "1e-6", ".2", "1", "1", "8"], log)
        command(["build/velocity_panels", "angular", "panels_probes", "1e-6", ".2", "1", ".01", "32"], log)
        command(["build/check_transport"], log)
        for name in [*METHODS, "base_original", "ref_original", "radial_truth", "panels_tight_radial", "ft07_refined_radial"]:
            command([sys.executable, "py/cmp_od.py", name], log)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--analyze", action="store_true")
    args = parser.parse_args()
    if not args.analyze:
        run()
    # Kept separate so analysis never recompiles or overwrites the truth table.
    command([sys.executable, "src/bench/write_summary.py"])
