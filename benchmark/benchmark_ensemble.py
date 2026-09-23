#!/usr/bin/env python3
"""Paired wall-clock timings of the complete executable; same seeds and inputs."""
import argparse
import csv
import statistics
import subprocess
import time
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--repeats', type=int, default=5)
parser.add_argument('--photons', type=int, default=128)
args = parser.parse_args()
if args.repeats < 1 or args.photons < 1:
    parser.error('positive repeats and photons required')
root = Path(__file__).resolve().parents[1]
rows = []
for beta0, mode in [(-.75,'O'),(-.75,'E'),(-.2,'O'),(-.2,'E')]:
    timing = {m: [] for m in ['event','fast']}
    outcomes = {}
    def run(method):
        begin = time.perf_counter()
        text = subprocess.check_output([str(root/'build/main'), '--method', method,
                 '--photons', str(args.photons), '--seed', '1234', '--beta0', str(beta0),
                 '--mode', mode], cwd=root, text=True)
        elapsed = time.perf_counter()-begin
        fields = dict(token.split('=',1) for token in text.split() if '=' in token)
        return elapsed, {k:int(fields[k]) for k in ['escaped','absorbed','scatterings']}
    for m in timing:
        run(m)
    # Reverse order on alternating repetitions to reduce timing-order bias.
    for repetition in range(args.repeats):
        for method in (['event','fast'] if repetition%2==0 else ['fast','event']):
            elapsed, outcome = run(method)
            timing[method].append(elapsed)
            outcomes[method] = outcome
            rows.append(dict(beta0=beta0, mode=mode, photons=args.photons,
                             repetition=repetition, method=method, seconds=elapsed, **outcome))
    old = statistics.median(timing['event']); fast = statistics.median(timing['fast'])
    print(f'beta0={beta0} mode={mode}: event={old:.6f}s fast={fast:.6f}s speedup={old/fast:.2f}x; '
          f'counts_match={outcomes["event"]==outcomes["fast"]}', flush=True)
    if outcomes['event'] != outcomes['fast']:
        print('Counts differ: investigate together with path-wise accuracy tests; counts alone are not an accuracy test.')
with (root/'benchmark/ensemble_speed.csv').open('w') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
