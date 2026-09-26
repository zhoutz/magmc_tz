"""Reproducible optical-depth comparison; run from the repository root.

Only generated copies in build/qagp_compare are adapted/instrumented.
The converged nonradial reference is extracted from the pinned local commit.
Requires g++-16 and Homebrew GSL, but no third-party Python packages.
"""
from pathlib import Path
import hashlib
import json
import math
import random
import shutil
import statistics
import subprocess

OUT = Path('build/qagp_compare')
REFERENCE = 'a907a14c98e03c8183a514791c2303a302d15bbe'
FLAGS = ['-std=c++23', '-O3', '-lgsl', '-I/opt/homebrew/include', '-L/opt/homebrew/lib']
HEADER = '''#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <limits>
#include <cstddef>
inline size_t bc_steps=0,bc_density=0,bc_quad=0,bc_geo=0;
inline bool bc_matched=false;
'''
MAIN = '''
int main(int argc,char**argv) {
 gsl_set_error_handler_off();
 std::ifstream in(argv[1]); std::ofstream out(argv[2]); out<<std::setprecision(17);
 int id,p; double b,m,a,z,e,r; std::string mode=argc>3?argv[3]:"default";
 while(in>>id>>b>>m>>a>>z>>e>>r>>p) {
  double tau=std::numeric_limits<double>::quiet_NaN(); int fail=0;
  std::string error;
  bc_steps=bc_density=bc_quad=bc_geo=0;
  RESET
  auto start=std::chrono::steady_clock::now();
  try { CALC } catch(std::exception const&ex) { fail=1; error=ex.what(); }
  double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  if(fail) std::cerr<<"case "<<id<<": "<<error<<'\\n';
  out<<id<<' '<<tau<<' '<<sec<<' '<<fail<<' '<<COUNTS<<'\\n';
 }
}
'''


def replace(text, old, new):
    assert old in text, old
    return text.replace(old, new)


def setup():
    OUT.mkdir(parents=True, exist_ok=True)
    for method in ['human', 'test_qagp', 'ref']:
        dest = OUT / method
        shutil.copytree('src', dest / 'src', dirs_exist_ok=True)
        source = dest / 'src/bench' / (method + '.cpp')
        s = source.read_text().split('int main()')[0]
        s = replace(s, 'Polarization pol, int n_knots',
                    'Polarization pol, double r0, double alpha, double az, int n_knots') if method != 'ref' else replace(s, 'Polarization pol) {\n  Boltzmann', 'Polarization pol, double r0, double alpha, double az) {\n  Boltzmann')
        s = replace(s, 'double3 n{0, 1, 0};', '''double3 theta{muz, 0, -r_hat.x}, phi{0, 1, 0};
  double3 e2 = std::cos(az)*theta + std::sin(az)*phi;
  double3 n = cross(r_hat, e2);''')
        s = replace(s, '.e2 = cross(n, r_hat),', '.e2 = e2,')
        s = replace(s, '{R_star, 0.0, 0.0}', '{r0, 0.0, alpha}')
        if method != 'ref':
            s = replace(s, 'Geometry geo(YVector const &y) const {', 'Geometry geo(YVector const &y) const {\n    ++bc_geo;')
            s = replace(s, 'auto integrand = [&](double path_length) {', 'auto integrand = [&](double path_length) {\n          ++bc_density;')
            s = replace(s, 'stepper.do_step(0.1 * stepper.y_old[0]);', 'stepper.do_step(0.1 * stepper.y_old[0]);\n    ++bc_steps;')
        else:
            s = replace(s, 'StepperDopr5<3, PhotonEvolution> stepper(photon_evolution, 1e-6, 1e-6);', 'StepperDopr5<3, PhotonEvolution> stepper(photon_evolution, bc_matched?1e-10:1e-6, bc_matched?1e-10:1e-6);')
            # Keep the original radial shortcut. Nonradial crossings can re-enter resonance.
            s = replace(s, 'bfield.calc_B(r, muz).length() * std::sqrt(1 - rs / r) / oi - 1;', 'bfield.calc_B(r, alpha == 0 ? muz : (std::cos(y[1])*r_hat + std::sin(y[1])*e2).z).length() * std::sqrt(1 - rs / r) / oi - 1;')
            s = replace(s, 'stepper.init(0.0,', 'if(alpha != 0) stepper.add_event([](double, YVector const& y){return y[0]-R_star;});\n  stepper.init(0.0,')
            s = replace(s, 'if (event_id == 0 || event_id == 1 || (event_id == 2 && stop_at_beta_zero)) break;', 'if (event_id == 0 || (alpha != 0 ? event_id == 3 : (event_id == 1 || (event_id == 2 && stop_at_beta_zero)))) break;')
            s = replace(s, '1e-6, 1e-6);', 'bc_matched?1e-8:1e-6, bc_matched?1e-8:1e-6);')
            s = replace(s, 'stepper.do_step(1e-2 * stepper.y_old[0]);', 'stepper.do_step(1e-2 * stepper.y_old[0]);\n    ++bc_steps;')
            s = replace(s, 'auto [r, psi, alpha] = stepper.dense_out(x);', '++bc_density;\n          auto [r, psi, alpha] = stepper.dense_out(x);')
        source.write_text(s)
        for header in ['quad.hpp', 'qags.hpp']:
            p = dest / 'src' / header
            t = p.read_text().replace('std::abort();', 'throw std::runtime_error("quadrature failure (see GSL status)");')
            t = t.replace('gsl_function gf;', '++bc_quad;\n    gsl_function gf;')
            t = t.replace('gsl_function integrand{};', '++bc_quad;\n  gsl_function integrand{};')
            p.write_text(t)
        calc = 'bc_matched=(mode=="matched"); tau=total_optical_depth(b,m,e,p?Polarization::E:Polarization::O,r,a,z);'
        driver = HEADER + f'#include "src/bench/{method}.cpp"\n' + MAIN.replace('RESET', '').replace('CALC', calc).replace('COUNTS', "bc_geo<<' '<<bc_density<<' '<<bc_steps<<' '<<bc_quad<<' '<<0")
        (dest / 'driver.cpp').write_text(driver)

    dest = OUT / 'truth'
    files = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', REFERENCE, 'src'], text=True).splitlines()
    for name in files:
        p = dest / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(subprocess.check_output(['git', 'show', REFERENCE + ':' + name]))
    calc = '''od::Ray ray(m,a,z,e,p?Polarization::E:Polarization::O,r);
 auto dist=od::thermal(b); od::Settings cfg;
 double scale=mode=="fine"?.03125:.125;
 cfg.cap=.01*scale; cfg.orbit_tol=1e-11*scale; cfg.tol=1e-8*scale;
 cfg.scan=1; cfg.transformed=true; cfg.reference=true;
 tau=od::integrate(ray,dist,cfg,c);'''
    (dest / 'driver.cpp').write_text(HEADER + '#include "src/bench/od_common.hpp"\n' + MAIN.replace('RESET', 'od::Counts c;').replace('CALC', calc).replace('COUNTS', "c.geometry<<' '<<c.density<<' '<<c.steps<<' '<<c.quadrature<<' '<<c.warnings"))
    for method in ['human', 'test_qagp', 'ref', 'truth']:
        subprocess.run(['g++-16', str(OUT / method / 'driver.cpp'), '-o', str(OUT / (method + '_run')), *FLAGS], check=True)

    with (OUT / 'radial.dat').open('w') as f:
        i = 0
        for k in range(1, 10):
            for m in range(10):
                for e in [.01, .1, 1, 10, 100]:
                    for p in [1, 0]:
                        f.write(f'{i} {-k/10} {m/10} 0 0 {e} 10 {p}\n')
                        i += 1
    rng = random.Random(20260927)
    with (OUT / 'random.dat').open('w') as f:
        for i in range(500):
            b = -rng.uniform(.1, .9)
            m = rng.uniform(-.999, .999)
            z = rng.uniform(0, 2*math.pi)
            e = 10**rng.uniform(-2, 2)
            r = 10 if i < 250 else 10**rng.uniform(math.log10(10.01), 2)
            a = math.acos(rng.uniform(0, 1) if i < 250 else rng.uniform(-1, 1))
            for p in [1, 0]:
                f.write(f'{2*i+1-p} {b:.17g} {m:.17g} {a:.17g} {z:.17g} {e:.17g} {r:.17g} {p}\n')
    paths = [*Path('src').rglob('*.hpp'), *Path('src/bench').glob('*.cpp'), Path('table/bfield_t10.txt'), Path('table/bench_od.txt')]
    manifest = dict(commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(), reference=REFERENCE, seed=20260927, flags=FLAGS,
                    compiler=subprocess.check_output(['g++-16', '--version'], text=True),
                    sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))


def rows(path):
    return [list(map(float, line.split())) for line in path.read_text().splitlines()]


def velocity():
    # Independent radial integration in beta: no orbit stepping, quadratic
    # resonance solver, or Boltzmann/Bessel implementation from the candidates.
    source = HEADER + '''
#include "human/src/bfield.hpp"
#include "human/src/constants.hpp"
#include "human/src/quad.hpp"
inline BField field("table/bfield_t10.txt",1e14,10);
double velocity_depth(double b0,double muz,double energy,int pol,double tol) {
 const double rs=1.4*schwarzschild_radius_of_sun_in_km, q=field.p+2;
 auto bv=field.calc_B(10,muz); double mu=bv.x/bv.length();
 double c=bv.length()*B_to_omega/energy;
 auto x=[&](double r){return c*std::pow(10/r,q)*std::sqrt(1-rs/r);};
 double G=x(10);
 if(G<=1) return 0;
 double lo=-(G*G-1)/(mu+G*std::sqrt(G*G+mu*mu-1));
 double s0=std::sqrt((1-b0)*(1+b0));
 double a=s0*(1+s0)/(b0*b0);
 auto fn=[&](double beta){
  double s=std::sqrt((1-beta)*(1+beta));
  double target=(1-beta*mu)/s;
  double left=10,right=10000;
  for(int j=0;j<60;++j){
   double mid=(left+right)/2;
   if(x(mid)>target)left=mid;else right=mid;
  }
  double r=(left+right)/2, L=std::sqrt(1-rs/r), Q=q-.5*rs/(r-rs);
  double mp=(mu-beta)/(1-beta*mu), P=pol?.5:.5*mp*mp;
  double gm1=beta*beta/(s*(1+s));
  // f(beta)/abs(mean_beta): the scaled Bessel normalization cancels exactly.
  return a*std::exp(-a*gm1)/(s*s*s)*P*(1-beta*mu)/(L*Q);
 };
 Quad quad;
 return (field.p+1)*pi*field.Bphi_over_Btheta(muz)*quad.qags(fn,lo,0,tol,tol);
}
'''
    source += MAIN.replace('RESET', '').replace('CALC', 'tau=velocity_depth(b,m,e,p,mode=="fine"?1e-12:1e-10);').replace('COUNTS', "0<<' '<<0<<' '<<0<<' '<<bc_quad<<' '<<0")
    (OUT / 'velocity.cpp').write_text(source)
    subprocess.run(['g++-16', str(OUT / 'velocity.cpp'), '-o', str(OUT / 'velocity_run'), *FLAGS], check=True)
    for level in ['coarse', 'fine']:
        with (OUT / f'velocity_{level}.log').open('w') as log:
            subprocess.run([str(OUT/'velocity_run'),str(OUT/'radial.dat'),str(OUT/f'velocity_{level}.dat'),level],stderr=log,check=True)
    coarse=rows(OUT/'velocity_coarse.dat'); fine=rows(OUT/'velocity_fine.dat')
    assert len(coarse)==len(fine)==900 and not any(r[3] for r in coarse+fine)
    stats=dict(convergence=errors(coarse,[r[1] for r in fine]), vs_table=errors(fine,[r[-1] for r in rows(Path('table/bench_od.txt'))]))
    (OUT/'velocity_metrics.json').write_text(json.dumps(stats,indent=2))
    print('velocity',json.dumps(stats),flush=True)


def run():
    # Serial and rotating method order; repeat 0 is a full warmup, then 5 timed runs.
    methods = ['human', 'test_qagp', 'ref', 'ref_matched']
    for suite in ['radial', 'random']:
        for rep in range(6):
            order = methods[rep % 4:] + methods[:rep % 4]
            for method in order:
                name = f'{suite}_{method}_{rep}'
                exe = 'ref' if method == 'ref_matched' else method
                with (OUT / (name + '.log')).open('w') as log:
                    subprocess.run([str(OUT / (exe + '_run')), str(OUT / (suite + '.dat')), str(OUT / (name + '.dat')), 'matched' if method == 'ref_matched' else 'default'], stderr=log, check=True)
                rr = rows(OUT / (name + '.dat'))
                assert len(rr) == (900 if suite == 'radial' else 1000)
                print(name, 'seconds', sum(r[2] for r in rr), 'failures', sum(r[3] for r in rr), flush=True)
        for level in ['coarse', 'fine']:
            name = f'{suite}_truth_{level}'
            with (OUT / (name + '.log')).open('w') as log:
                subprocess.run([str(OUT / 'truth_run'), str(OUT / (suite + '.dat')), str(OUT / (name + '.dat')), level], stderr=log, check=True)
            rr = rows(OUT / (name + '.dat'))
            print(name, 'seconds', sum(r[2] for r in rr), 'failures', sum(r[3] for r in rr), flush=True)


def errors(actual, truth, ids=None):
    ids = range(len(actual)) if ids is None else ids
    pairs = [(i, abs(actual[i][1]-truth[i]), abs(actual[i][1]-truth[i])/abs(truth[i]) if abs(truth[i]) > 1e-10 else None)
             for i in ids if not actual[i][3] and math.isfinite(truth[i])]
    relative = sorted(x[2] for x in pairs if x[2] is not None)
    return dict(n=len(pairs), max_abs=max(x[1] for x in pairs), max_rel=max(relative),
                p99_rel=relative[math.ceil(.99*len(relative))-1], median_rel=statistics.median(relative),
                worst_rel=max((x for x in pairs if x[2] is not None), key=lambda x:x[2]),
                above_1e6=sum(v > 1e-6 for v in relative), above_1e8=sum(v > 1e-8 for v in relative))


def analyze():
    result = {}
    for suite in ['radial', 'random']:
        coarse = rows(OUT / f'{suite}_truth_coarse.dat')
        fine = rows(OUT / f'{suite}_truth_fine.dat')
        assert len(coarse) == len(fine) == (900 if suite == 'radial' else 1000)
        assert not any(r[3] for r in fine + coarse), 'Reference failures need investigation'
        truth = [r[1] for r in fine] if suite == 'random' else [r[-1] for r in rows(Path('table/bench_od.txt'))]
        result[suite] = {'convergence': errors(coarse, [r[1] for r in fine]), 'truth_vs_target': errors(fine, truth)}
        methods = ['human', 'test_qagp', 'ref', 'ref_matched']
        data = {m: [rows(OUT / f'{suite}_{m}_{rep}.dat') for rep in range(1, 6)] for m in methods}
        common = [i for i in range(len(truth)) if all(not data[m][0][i][3] for m in methods)]
        for m, runs in data.items():
            first = runs[0]
            for other in runs[1:]:
                assert all(a[3] == b[3] and (a[3] or a[1] == b[1]) for a, b in zip(first, other)), 'Nondeterministic outputs'
            times = [sum(r[2] for r in run) for run in runs]
            result[suite][m] = dict(errors=errors(first, truth), failures=sum(r[3] for r in first),
                seconds=statistics.median(times), seconds_min=min(times), seconds_max=max(times),
                common_errors=errors(first, truth, common), common_seconds=statistics.median(sum(run[i][2] for i in common) for run in runs),
                counts=dict(zip(['geometry','density','steps','quadrature','warnings'], [sum(r[j] for r in first) for j in range(4, 9)])))
        pair = [i for i in range(len(truth)) if not data['human'][0][i][3] and not data['test_qagp'][0][i][3]]
        result[suite]['human_qagp_common'] = {
            m: dict(errors=errors(data[m][0], truth, pair), seconds=statistics.median(sum(run[i][2] for i in pair) for run in data[m]))
            for m in ['human', 'test_qagp']}
        result[suite]['human_qagp_max_difference'] = max(abs(data['human'][0][i][1]-data['test_qagp'][0][i][1]) for i in pair)
        inputs = rows(OUT / (suite + '.dat'))
        for m in methods:
            result[suite][m]['failure_inputs'] = [inputs[i] for i in range(len(truth)) if data[m][0][i][3]]
            groups = {'E': [i for i, r in enumerate(inputs) if r[7] == 1], 'O': [i for i, r in enumerate(inputs) if r[7] == 0]}
            if suite == 'random':
                groups.update(surface=list(range(500)), magnetosphere=list(range(500,1000)))
            result[suite][m]['groups'] = {name: dict(errors=errors(data[m][0],truth,ids), failures=sum(data[m][0][i][3] for i in ids)) for name,ids in groups.items()}
    (OUT / 'metrics.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


def refine():
    coarse=rows(OUT/'random_truth_coarse.dat'); fine=rows(OUT/'random_truth_fine.dat')
    ids=set(sorted(range(1000),key=lambda i:abs(coarse[i][1]-fine[i][1]),reverse=True)[:12])
    for method in ['human','test_qagp','ref','ref_matched']:
        actual=rows(OUT/f'random_{method}_1.dat')
        good=[i for i in range(1000) if not actual[i][3]]
        ids.update(i for i in range(1000) if actual[i][3])
        ids.update(sorted(good,key=lambda i:abs(actual[i][1]-fine[i][1]),reverse=True)[:5])
        ids.update(sorted(good,key=lambda i:abs(actual[i][1]-fine[i][1])/max(abs(fine[i][1]),1e-10),reverse=True)[:5])
    inputs=(OUT/'random.dat').read_text().splitlines()
    (OUT/'refinement.dat').write_text(''.join(inputs[i]+'\n' for i in sorted(ids)))
    driver=(OUT/'truth/driver.cpp').read_text()
    driver=replace(driver,'mode=="fine"?.03125:.125','.0078125')
    driver=replace(driver,'cfg.tol=1e-8*scale;','cfg.tol=1e-12;')
    common=OUT/'truth/src/bench/od_common.hpp'
    # At 1/4 of the fine step cap, a complete trajectory can exceed 100k steps.
    (common.parent/'od_common_refinement.hpp').write_text(common.read_text().replace('it < 100000','it < 500000'))
    driver=replace(driver,'src/bench/od_common.hpp','src/bench/od_common_refinement.hpp')
    (OUT/'truth/refinement.cpp').write_text(driver)
    subprocess.run(['g++-16',str(OUT/'truth/refinement.cpp'),'-o',str(OUT/'refinement_run'),*FLAGS],check=True)
    with (OUT/'refinement.log').open('w') as log:
        subprocess.run([str(OUT/'refinement_run'),str(OUT/'refinement.dat'),str(OUT/'refinement_result.dat')],stderr=log,check=True)
    rr=rows(OUT/'refinement_result.dat')
    assert not any(r[3] for r in rr)
    summary=[dict(id=int(r[0]),tau=r[1],abs_change=abs(r[1]-fine[int(r[0])][1]),rel_change=abs(r[1]-fine[int(r[0])][1])/max(abs(r[1]),1e-10),warnings=r[8]) for r in rr]
    (OUT/'refinement_metrics.json').write_text(json.dumps(summary,indent=2))
    print('refinement',len(rr),'max_abs',max(r['abs_change'] for r in summary),'max_rel',max(r['rel_change'] for r in summary))


def originals():
    # Verify that the general initial conditions and failure instrumentation
    # have not changed any of the original radial results.
    dest=OUT/'pristine'
    shutil.copytree('src',dest/'src',dirs_exist_ok=True)
    (dest/'table').mkdir(exist_ok=True)
    (dest/'output').mkdir(exist_ok=True)
    shutil.copyfile('table/bfield_t10.txt',dest/'table/bfield_t10.txt')
    result={}
    for method in ['human','test_qagp','ref']:
        exe=(dest/method).resolve()
        subprocess.run(['g++-16',str(dest/'src/bench'/f'{method}.cpp'),'-o',str(exe),*FLAGS],check=True)
        subprocess.run([str(exe)],cwd=dest,check=True)
        raw=rows(dest/'output'/f'{method}.txt')
        adapted=rows(OUT/f'radial_{method}_1.dat')
        assert len(raw)==len(adapted)==900
        result[method]=max(abs(a[-1]-b[1]) for a,b in zip(raw,adapted))
        assert result[method]<1e-12
    (OUT/'original_validation.json').write_text(json.dumps(result,indent=2))
    print('originals max absolute differences',result)


if __name__ == '__main__':
    import sys
    stages = sys.argv[1:] or ['setup', 'run', 'velocity', 'analyze', 'refine', 'originals']
    for stage in stages:
        globals()[stage]()
