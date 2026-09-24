#!/usr/bin/env python3
"""Reproduce all general-angle experiments. No third-party Python packages."""
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

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output'
NAMES=['baseline','ft07','phase_cap','event_guard']
RUNS=[('baseline','1e-6','.1','.01'),('ft07','1e-6','.01','.01'),
      ('phase_cap','1e-6','.1','.01'),('event_guard','1e-6','.1','.01'),
      ('baseline','1e-8','.1','.01'),('ft07','1e-8','.01','.01'),
      ('phase_cap','1e-8','.1','.01'),('event_guard','1e-8','.1','.01'),
      ('phase_cap','1e-6','.03','.01'),('phase_cap','1e-6','.3','.01'),
      ('event_guard','1e-6','.025','.01'),('event_guard','1e-6','.05','.01')]
RUNS += [(n,'1e-6','.01' if n=='ft07' else '.1',h) for n in ['baseline','phase_cap','event_guard'] for h in ['.001','1']]

def filename(run):
    name,tol,res,h=run
    return OUT/f'{name}_tol{tol}_p{res}_h{h}.csv'

def summarize():
    summary=[]
    for run in RUNS:
        path=filename(run)
        if not path.exists():continue
        with path.open() as f:
            meta=dict(w.split('=',1) for w in f.readline().removeprefix('# ').split())
            rows=list(csv.DictReader(f))
        good=[r for r in rows if r['status']=='ok'];errors=sorted(float(r['rel_error']) for r in good)
        summary.append(dict(file=path.name,**meta,cases=len(rows),failures=len(rows)-len(good),
             maximum_relative=max(errors),median_relative=statistics.median(errors),p95_relative=errors[math.ceil(.95*len(errors))-1],
             maximum_absolute=max(float(r['abs_error']) for r in good),above_1e3=sum(e>1e-3 for e in errors),
             missed_90pct=sum(float(r['tau'])<.1*float(r['reference']) for r in good),
             mean_evaluations=statistics.mean(int(r['evaluations']) for r in good),
             mean_steps=statistics.mean(int(r['steps']) for r in good)))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (OUT/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=summary[0]);w.writeheader();w.writerows(summary)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repeats',type=int,default=5)
    parser.add_argument('--compiler',default=os.environ.get('CXX',shutil.which('g++-16') or 'c++'))
    parser.add_argument('--methods',nargs='+',choices=NAMES);parser.add_argument('--skip-validation',action='store_true')
    args=parser.parse_args();os.chdir(ROOT)
    flags=['-std=c++23','-O3','-Wall','-Wextra','-Wno-psabi','-pedantic']
    selected=args.methods or NAMES
    metadata=dict(date=time.strftime('%Y-%m-%d %H:%M:%S %z'),platform=platform.platform(),compiler=subprocess.check_output([args.compiler,'--version'],text=True),flags=flags,repeats=args.repeats,
        reference_sha256=hashlib.sha256(Path('table/bench_od.txt').read_bytes()).hexdigest(),field_sha256=hashlib.sha256(Path('table/bfield_t10.txt').read_bytes()).hexdigest())
    (OUT/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    with (OUT/'run.log').open('a' if args.methods else 'w') as log:
        def run(command):
            log.write('$ '+' '.join(map(str,command))+'\n');log.flush()
            p=subprocess.run(command,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            log.write(p.stdout+'\n');log.flush();print(p.stdout,end='',flush=True);p.check_returncode()
        for name in selected:run([args.compiler,*flags,f'benchmark/{name}.cpp','-o',f'output/{name}'])
        for config in RUNS:
            name,tol,res,h=config
            if name not in selected:continue
            run([f'output/{name}',str(filename(config)),tol,res,str(args.repeats),h])
            summarize()
        if not args.skip_validation:
            for name in ['angular_validation','transport_validation']:
                run([args.compiler,*flags,f'benchmark/{name}.cpp','-o',f'output/{name}'])
                run([f'output/{name}'])
        run([args.compiler,*flags,'src/main.cpp','-o','output/main'])
        with (OUT/'main_results.txt').open('w') as f:subprocess.run(['output/main'],stdout=f,check=True)
    summarize()

if __name__=='__main__':main()
