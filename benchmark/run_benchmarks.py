#!/usr/bin/env python3
"""Run from any cwd. Standard library only; all generated files stay in output/."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"

def run(command, log):
    p = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.write("$ " + " ".join(map(str, command)) + "\n" + p.stdout + "\n")
    log.flush()
    print(p.stdout, end="", flush=True)
    p.check_returncode()

def summarize(paths):
    summary=[]
    for path in paths:
        with path.open() as f:
            metadata=dict(x.split("=",1) for x in f.readline().removeprefix("# ").split())
            rows=list(csv.DictReader(f))
        good=[r for r in rows if r['status']=='ok']
        errors=sorted(float(r['rel_error']) for r in good)
        worst=max(good,key=lambda r:float(r['rel_error']))
        summary.append(dict(file=path.name, **metadata, failures=len(rows)-len(good),
            max_relative=max(errors), median_relative=statistics.median(errors), p95_relative=errors[math.ceil(.95*len(errors))-1],
            max_absolute=max(float(r['abs_error']) for r in good),
            above_1e3=sum(e>1e-3 for e in errors), above_1e4=sum(e>1e-4 for e in errors),
            missed_90pct=sum(float(r['tau'])<.1*float(r['reference']) for r in good),
            mean_evaluations=statistics.mean(int(r['evaluations']) for r in good),
            mean_steps=statistics.mean(int(r['steps']) for r in good),
            mean_rejections=statistics.mean(int(r['rejected']) for r in good),
            worst_case="/".join(worst[k] for k in ['b0','muz','pol','omega_inf'])))
    with (OUT/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=summary[0].keys());w.writeheader();w.writerows(summary)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repeats',type=int,default=7)
    parser.add_argument('--compiler',default=os.environ.get('CXX',shutil.which('g++-16') or 'c++'))
    parser.add_argument('--methods',nargs='+',help='Re-run only these methods; keep other existing CSV measurements')
    args=parser.parse_args()
    os.chdir(ROOT)
    names=['ft07','baseline','thermal_cap','spatial_regularized','velocity']
    flags=['-std=c++23','-O3','-Wall','-Wextra','-pedantic']
    metadata=dict(date=time.strftime('%Y-%m-%d %H:%M:%S %z'),platform=platform.platform(),machine=platform.machine(),
                  compiler=subprocess.check_output([args.compiler,'--version'],text=True),flags=flags,
                  reference_sha256=hashlib.sha256((ROOT/'table/bench_od.txt').read_bytes()).hexdigest(),
                  field_sha256=hashlib.sha256((ROOT/'table/bfield_t10.txt').read_bytes()).hexdigest(),repeats=args.repeats)
    (OUT/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    selected=args.methods or names
    if any(name not in names for name in selected):parser.error('Unknown method')
    with (OUT/'run.log').open('a' if args.methods else 'w') as log:
        for name in selected+['validation']:
            run([args.compiler,*flags,'benchmark/'+name+'.cpp','-o','output/'+name],log)
        runs=[('ft07','1e-6','.01','.01'),('baseline','1e-6','.01','.01'),('thermal_cap','1e-6','.01','.01'),
              ('spatial_regularized','1e-6','.01','.01'),('velocity','1e-6','.01','.01')]
        runs += [(name,tol,'.01','.01') for name in names for tol in ['1e-7','1e-8','1e-10']]
        runs += [('ft07','1e-8',p,'.01') for p in ['.005','.02']]
        runs += [('thermal_cap','1e-7',p,'.01') for p in ['.002','.05','.1']]
        # Initial-step sensitivity: same physics and table, different sample phases.
        runs += [(name,'1e-6','.01',h) for name in ['baseline','thermal_cap'] for h in ['.001','.1','1']]
        paths=[]
        for name,tol,parameter,h in runs:
            path=OUT/f'{name}_tol{tol}_p{parameter}_h{h}.csv';paths.append(path)
            if name not in selected:
                if not path.exists():raise RuntimeError('Missing measurement: '+str(path))
                continue
            run(['output/'+name,str(path),tol,parameter,str(args.repeats),h],log)
        run(['output/validation'],log)
        run([args.compiler,*flags,'src/main.cpp','-o','output/main'],log)
        with (OUT/'main_results.txt').open('w') as f: subprocess.run(['output/main'],stdout=f,check=True)
    summarize(paths)

if __name__=='__main__':main()
