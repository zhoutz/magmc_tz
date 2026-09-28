"""Compare current human.cpp and human.1.cpp using the established 1900 rays.

Run from the repository root. Numerical sources are changed only in build copies.
Requires reference artifacts from compare_all_optical_depth.py.
"""
from pathlib import Path
import hashlib
import json
import shutil
import statistics
import subprocess
import sys

import compare_qagp as base

OUT = Path('build/human1_compare')
PREVIOUS = Path('build/optical_compare')
METHODS = ['human', 'human.1']


def setup():
    OUT.mkdir(parents=True, exist_ok=True)
    # The reused reference has the same field table and physical/numerical helpers.
    prior = json.loads((PREVIOUS/'manifest.json').read_text())['sha256']
    for name in ['table/bfield_t10.txt', 'table/bench_od.txt', 'src/constants.hpp',
                 'src/bessel_K1_scaled.hpp', 'src/doubles.hpp', 'src/hunt.hpp', 'src/dopr5.hpp', 'src/roots.hpp']:
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == prior[name], name
    for suite in ['radial', 'random']:
        shutil.copyfile(PREVIOUS/(suite+'.dat'), OUT/(suite+'.dat'))
    for method in METHODS:
        dest=OUT/method
        shutil.copytree('src',dest/'src',dirs_exist_ok=True)
        path=dest/'src/bench'/f'{method}.cpp'
        s=path.read_text().split('int main()')[0]
        s=base.replace(s,'Polarization pol, int n_knots','Polarization pol, double r0, double alpha, double az, int n_knots')
        s=base.replace(s,'double3 n{0, 1, 0};','''double3 theta{muz, 0, -r_hat.x}, phi{0, 1, 0};
  double3 e2=std::cos(az)*theta+std::sin(az)*phi;
  double3 n=cross(r_hat,e2);''')
        s=base.replace(s,'.e2 = cross(n, r_hat),','.e2 = e2,')
        s=base.replace(s,'{R_star, 0.0, 0.0}','{r0, 0.0, alpha}')
        s=base.replace(s,'Geometry geo(YVector const &y) const {','Geometry geo(YVector const &y) const {\n    ++bc_geo;')
        s=base.replace(s,'auto integrand = [&](double path_length) {','auto integrand = [&](double path_length) {\n    ++bc_density;')
        s=base.replace(s,'stepper.do_step(0.1 * stepper.y_old[0]);','stepper.do_step(0.1 * stepper.y_old[0]);\n    ++bc_steps;')
        if method=='human.1':
            s=base.replace(s,'if (fb.f(betas[0]) == 0 && fb.f(betas[1]) == 0) continue;', '''if (fb.f(betas[0]) == 0 && fb.f(betas[1]) == 0) {
          bool in_support=false;
          for(double beta:betas) in_support |= beta>fb.b_min && beta<fb.b_max;
          if(in_support) ++bc_underflow; else ++bc_support;
          continue;
        }''')
        path.write_text(s)
        q=dest/'src/quad.hpp'
        q.write_text(q.read_text().replace('gsl_function gf;','++bc_quad;\n    gsl_function gf;').replace('std::abort();','throw std::runtime_error("quadrature failure");'))
        driver=base.HEADER+'inline size_t bc_support=0,bc_underflow=0;\n'+f'#include "src/bench/{method}.cpp"\n'
        driver+=base.MAIN.replace('RESET','bc_support=bc_underflow=0;').replace('CALC','tau=total_optical_depth(b,m,e,p?Polarization::E:Polarization::O,r,a,z);').replace('COUNTS',"bc_geo<<' '<<bc_density<<' '<<bc_steps<<' '<<bc_quad<<' '<<0<<' '<<bc_support<<' '<<bc_underflow")
        (dest/'driver.cpp').write_text(driver)
        subprocess.run(['g++-16',str(dest/'driver.cpp'),'-o',str(OUT/(method+'_run')),*base.FLAGS],check=True)
    paths=[*Path('src').rglob('*.hpp'), *(Path('src/bench')/(m+'.cpp') for m in METHODS), Path('table/bfield_t10.txt'), Path('table/bench_od.txt'), PREVIOUS/'random_truth_fine.dat']
    manifest=dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        flags=base.FLAGS, compiler=subprocess.check_output(['g++-16','--version'],text=True),seed=20260927,
        sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))


def run():
    for suite in ['radial','random']:
        # One warmup, then seven measurements with alternating execution order.
        for rep in range(8):
            for m in METHODS[::1 if rep%2==0 else -1]:
                stem=f'{suite}_{m}_{rep}'
                with (OUT/(stem+'.log')).open('w') as log:
                    subprocess.run([str(OUT/(m+'_run')),str(OUT/(suite+'.dat')),str(OUT/(stem+'.dat'))],stderr=log,check=True)
                rr=base.rows(OUT/(stem+'.dat'))
                assert len(rr)==(900 if suite=='radial' else 1000)
                print(stem,'seconds',sum(r[2] for r in rr),'failures',sum(r[3] for r in rr),flush=True)


def analyze():
    result={}
    for suite in ['radial','random']:
        truth=[r[-1] for r in base.rows(Path('table/bench_od.txt'))] if suite=='radial' else [r[1] for r in base.rows(PREVIOUS/'random_truth_fine.dat')]
        batches={m:[base.rows(OUT/f'{suite}_{m}_{rep}.dat') for rep in range(1,8)] for m in METHODS}
        common=[i for i in range(len(truth)) if all(not batches[m][0][i][3] for m in METHODS)]
        inputs=base.rows(OUT/(suite+'.dat'))
        d={}
        for m,runs in batches.items():
            first=runs[0]
            for other in runs[1:]:
                assert all(a[3]==b[3] and (a[3] or a[1]==b[1]) for a,b in zip(first,other))
            times=[sum(r[2] for r in run) for run in runs]
            groups={'E':[i for i,r in enumerate(inputs) if r[-1]==1], 'O':[i for i,r in enumerate(inputs) if r[-1]==0]}
            if suite=='random':groups.update(surface=list(range(500)),magnetosphere=list(range(500,1000)))
            d[m]=dict(seconds=statistics.median(times),min_seconds=min(times),max_seconds=max(times),
                common_seconds=statistics.median(sum(run[i][2] for i in common) for run in runs),
                errors=base.errors(first,truth),common_errors=base.errors(first,truth,common),
                failures=sum(r[3] for r in first),failure_ids=[int(r[0]) for r in first if r[3]],
                counts=dict(zip(['geometry','density','steps','quadrature','warnings','support_skips','underflow_skips'],[sum(r[k] for r in first) for k in range(4,11)])),
                groups={name:dict(errors=base.errors(first,truth,ids),seconds=statistics.median(sum(run[i][2] for i in ids) for run in runs)) for name,ids in groups.items()})
        d['speedup']=d['human']['seconds']/d['human.1']['seconds']
        d['time_reduction_percent']=100*(1-d['human.1']['seconds']/d['human']['seconds'])
        d['common_speedup']=d['human']['common_seconds']/d['human.1']['common_seconds']
        d['max_difference']=max(abs(batches['human'][0][i][1]-batches['human.1'][0][i][1]) for i in common)
        result[suite]=d
    (OUT/'metrics.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


def originals():
    result={}
    for m in METHODS:
        dest=OUT/'pristine'/m
        shutil.copytree('src',dest/'src',dirs_exist_ok=True)
        (dest/'output').mkdir(exist_ok=True);(dest/'table').mkdir(exist_ok=True)
        shutil.copyfile('table/bfield_t10.txt',dest/'table/bfield_t10.txt')
        exe=(dest/'run').resolve()
        subprocess.run(['g++-16',str(dest/'src/bench'/f'{m}.cpp'),'-o',str(exe),*base.FLAGS],check=True)
        subprocess.run([str(exe)],cwd=dest,check=True)
        rr=base.rows(dest/'output/human.txt'); adapted=base.rows(OUT/f'radial_{m}_1.dat')
        assert len(rr)==len(adapted)==900
        result[m]=max(abs(r[-1]-a[1]) for r,a in zip(rr,adapted))
        assert result[m]==0
    (OUT/'original_validation.json').write_text(json.dumps(result,indent=2))


if __name__=='__main__':
    for stage in sys.argv[1:] or ['setup','run','analyze','originals']:
        globals()[stage]()
