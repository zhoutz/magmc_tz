"""Compare six historical optical-depth methods on identical inputs.

Run from the repository root. Preserves the earlier build/qagp_compare results.
"""
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys

import compare_qagp as base

OUT = Path('build/optical_compare')
base.OUT = OUT
NEW = ['global_sqrt', 'velocity_panels', 'transport_fast']
PRIMARY = ['human', 'ref', 'test_qagp', *NEW]
CONFIGS = {
    **{m: (m, 'default') for m in PRIMARY},
    'ref_matched': ('ref', 'matched'),
    'velocity_panels_matched': ('velocity_panels', 'matched'),
    'transport_fast_matched': ('transport_fast', 'matched'),
    'velocity_panels_cached': ('velocity_panels', 'cached'),
}
COMMITS = {'global_sqrt': 'a907a14c98e03c8183a514791c2303a302d15bbe',
           'velocity_panels': '9204013bc5919c2edad4e53c3c58bfe8cde0a4db',
           'transport_fast': '4be9c249402dab313c171bc04244bc2be6a39563'}


def setup():
    base.setup()
    for m in NEW:
        subprocess.run(['g++-16', f'src/bench/{m}.cpp', '-o', str(OUT / (m + '_run')), *base.FLAGS], check=True)
    manifest = json.loads((OUT / 'manifest.json').read_text())
    manifest.update(ports=COMMITS, configurations=CONFIGS)
    for name in ['py/compare_all_optical_depth.py', 'py/compare_qagp.py']:
        manifest['sha256'][name] = hashlib.sha256(Path(name).read_bytes()).hexdigest()
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))


def execute(suite, method, rep, config=None):
    exe, mode = config or CONFIGS[method]
    stem = f'{suite}_{method}_{rep}'
    with (OUT / (stem + '.log')).open('w') as log:
        subprocess.run([str(OUT / (exe + '_run')), str(OUT / (suite + '.dat')),
                        str(OUT / (stem + '.dat')), mode], stderr=log, check=True)
    rr = base.rows(OUT / (stem + '.dat'))
    assert len(rr) == (900 if suite == 'radial' else 1000)
    print(stem, 'seconds', sum(r[2] for r in rr), 'failures', sum(r[3] for r in rr), flush=True)


def run():
    methods = list(CONFIGS)
    for suite in ['radial', 'random']:
        for rep in range(6):
            k = rep % len(methods)
            for method in methods[k:] + methods[:k]:
                execute(suite, method, rep)
        for level in ['coarse', 'fine']:
            execute(suite, 'truth', level, ('truth', level))


def analyze():
    results = {}
    for suite in ['radial', 'random']:
        fine = base.rows(OUT / f'{suite}_truth_fine.dat')
        coarse = base.rows(OUT / f'{suite}_truth_coarse.dat')
        assert not any(r[3] for r in coarse + fine)
        truth = [r[1] for r in fine] if suite == 'random' else [r[-1] for r in base.rows(Path('table/bench_od.txt'))]
        inputs = base.rows(OUT / (suite + '.dat'))
        runs = {m: [base.rows(OUT / f'{suite}_{m}_{rep}.dat') for rep in range(1, 6)] for m in CONFIGS}
        common = [i for i in range(len(truth)) if all(not runs[m][0][i][3] for m in PRIMARY)]
        result = {'reference_convergence': base.errors(coarse, [r[1] for r in fine]), 'common_ids': common}
        for method, batches in runs.items():
            rr = batches[0]
            for other in batches[1:]:
                assert all(a[3] == b[3] and (a[3] or a[1] == b[1]) for a, b in zip(rr, other)), 'Nondeterministic result'
            times = [sum(r[2] for r in b) for b in batches]
            groups = {'E': [i for i, r in enumerate(inputs) if r[7] == 1],
                      'O': [i for i, r in enumerate(inputs) if r[7] == 0]}
            if suite == 'random':
                groups.update(surface=list(range(500)), magnetosphere=list(range(500, 1000)))
            setup = []
            for rep in range(1, 6):
                setup.append(sum(float(line.split('=')[1]) for line in (OUT / f'{suite}_{method}_{rep}.log').read_text().splitlines()
                                 if line.startswith('preparation_seconds=')))
            result[method] = dict(seconds=statistics.median(times), min_seconds=min(times), max_seconds=max(times),
                failures=sum(r[3] for r in rr), errors=base.errors(rr, truth),
                common_errors=base.errors(rr, truth, common),
                common_seconds=statistics.median(sum(b[i][2] for i in common) for b in batches),
                failure_inputs=[inputs[i] for i, r in enumerate(rr) if r[3]],
                groups={g: dict(errors=base.errors(rr, truth, ids), failures=sum(rr[i][3] for i in ids)) for g, ids in groups.items()},
                counts=dict(zip(['geometry', 'density', 'steps', 'quadrature', 'warnings'], [sum(r[k] for r in rr) for k in range(4, 9)])),
                timed_preparation_seconds=statistics.median(sum(r[9] if len(r) > 9 else 0 for r in b) for b in batches),
                cached_preparation_seconds=statistics.median(setup),
                missed_nonzero=sum(not r[3] and r[1] == 0 and abs(truth[i]) > 1e-10 for i, r in enumerate(rr)))
        results[suite] = result
    (OUT / 'metrics.json').write_text(json.dumps(results, indent=2))
    for suite, result in results.items():
        print(suite, 'common', len(result['common_ids']))
        for m in CONFIGS:
            r = result[m]
            print(m, r['seconds'], 'fail', r['failures'], 'max_abs', r['errors']['max_abs'],
                  'max_rel', r['errors']['max_rel'], 'p99', r['errors']['p99_rel'])


def refine():
    base.refine(list(CONFIGS))


def audit():
    # Supplemental convergence/configuration checks, outside the timing ranking.
    for suite in ['radial', 'random']:
        for method in NEW:
            execute(suite, method + '_audit', 0, (method, 'audit'))
        execute(suite, 'transport_fast_guarded', 0, ('transport_fast', 'guarded'))


def audit_analyze():
    result = {}
    for suite in ['radial', 'random']:
        truth = [r[1] for r in base.rows(OUT/f'{suite}_truth_fine.dat')] if suite == 'random' else [r[-1] for r in base.rows(Path('table/bench_od.txt'))]
        result[suite] = {}
        for method in [m+'_audit' for m in NEW] + ['transport_fast_guarded']:
            rr=base.rows(OUT/f'{suite}_{method}_0.dat')
            result[suite][method]=dict(errors=base.errors(rr,truth), failures=sum(r[3] for r in rr),
                seconds=sum(r[2] for r in rr), failure_ids=[int(r[0]) for r in rr if r[3]],
                warnings=sum(r[8] for r in rr))
    (OUT/'audit_metrics.json').write_text(json.dumps(result,indent=2))
    print('audit metrics',json.dumps(result),flush=True)


def validate_ports():
    """Run unmodified historical numerical code with the same 1900 inputs."""
    snippets = {
        'global_sqrt': ('#include "src/bench/od_global.hpp"\n', '''
 od::Ray ray(m,a,z,e,p?Polarization::E:Polarization::O,r);
 auto dist=od::thermal(b);od::Settings cfg;od::Counts c;
 cfg.cap=.5;cfg.orbit_tol=1e-10;cfg.tol=1e-8;cfg.scan=1;
 tau=od::global_integrate(ray,dist,cfg,c);'''),
        'velocity_panels': ('#include "src/bench/bench_transport.hpp"\ninline BField field("table/bfield_t10.txt",1e14,10);\n', '''
 auto dist=transport::boltzmann_distribution(b);
 auto photon=transport::make_photon(r,m,a,z,e,p?Polarization::E:Polarization::O);
 tau=transport::integrate(field,dist,photon,transport::Method::panels).tau;'''),
        'transport_fast': ('#include "src/transport_fast.hpp"\ninline BField field("table/bfield_t10.txt",1e14,10);\n', '''
 Boltzmann dist(b);double3 er{std::sqrt(1-m*m),0,m},theta{m,0,-er.x},phi{0,1,0};
 double3 tangent=std::cos(z)*theta+std::sin(z)*phi;
 Photon photon{cross(er,tangent),er,tangent,r,0,a,e,p?Polarization::E:Polarization::O};
 tau=fast_transport::propagate(field,dist,photon).tau;'''),
    }
    report = {}
    for method in NEW:
        dest = OUT / 'historical' / method
        commit = COMMITS[method]
        paths = subprocess.check_output(['git','ls-tree','-r','--name-only',commit,'src'],text=True).splitlines()
        for name in paths:
            target = dest / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(subprocess.check_output(['git','show',commit+':'+name]))
        include, calc = snippets[method]
        driver = base.HEADER + '#include <format>\n#include <gsl/gsl_errno.h>\n' + include + base.MAIN.replace('RESET','').replace('CALC',calc).replace('COUNTS',"0<<' '<<0<<' '<<0<<' '<<0<<' '<<0")
        (dest/'driver.cpp').write_text(driver)
        executable = dest/'run'
        subprocess.run(['g++-16',str(dest/'driver.cpp'),'-o',str(executable),*base.FLAGS],check=True)
        report[method] = {}
        for suite in ['radial','random']:
            target = dest/(suite+'.dat')
            with (dest/(suite+'.log')).open('w') as log:
                subprocess.run([str(executable),str(OUT/(suite+'.dat')),str(target)],stderr=log,check=True)
            old=base.rows(target);new=base.rows(OUT/f'{suite}_{method}_1.dat')
            assert len(old)==len(new)==(900 if suite=='radial' else 1000)
            assert all(a[3]==b[3] for a,b in zip(old,new)), 'Port changed failure behavior'
            differences=[abs(a[1]-b[1]) for a,b in zip(old,new) if not a[3]]
            relative=[abs(a[1]-b[1])/max(abs(a[1]),1e-10) for a,b in zip(old,new) if not a[3]]
            report[method][suite]=dict(max_abs=max(differences),max_rel=max(relative),failures=sum(r[3] for r in old))
            assert max(differences)<1e-10, 'Port changed optical depths'
        print('historical validation',method,report[method],flush=True)
    (OUT/'port_validation.json').write_text(json.dumps(report,indent=2))


def standalone():
    report = {}
    for method in NEW:
        subprocess.run([str(OUT/(method+'_run'))],check=True)
        actual=base.rows(Path('output')/(method+'.txt'))
        expected=base.rows(OUT/f'radial_{method}_1.dat')
        assert len(actual)==len(expected)==900
        report[method]=max(abs(a[-1]-b[1]) for a,b in zip(actual,expected))
        assert report[method]==0
    (OUT/'standalone_validation.json').write_text(json.dumps(report,indent=2))


def diagnose_warnings():
    rr=base.rows(OUT/'random_global_sqrt_1.dat')
    ids=[i for i,r in enumerate(rr) if r[8]]
    inputs=(OUT/'random.dat').read_text().splitlines()
    (OUT/'global_warning_inputs.dat').write_text(''.join(inputs[i]+'\n' for i in ids))
    if not ids:
        return
    source=Path('src/bench/global_sqrt.cpp').read_text()
    source=base.replace(source,'    if (status != GSL_SUCCESS) {',
                        '    int initial_status = status;\n    if (status != GSL_SUCCESS) {')
    marker='    if (!std::isfinite(value)) throw std::runtime_error("nonfinite quadrature");'
    source=base.replace(source,marker,
        '    if (initial_status != GSL_SUCCESS) std::fprintf(stderr, "initial_status=%d final_status=%d value=%.17g error=%.17g tol=%.17g\\n", initial_status,status,value,error,tol);\n'+marker)
    source=base.replace(source,'  od::Ray ray(c.muz','  std::cerr << "case " << c.id << "\\n";\n  od::Ray ray(c.muz')
    (OUT/'global_diagnostics.cpp').write_text(source)
    subprocess.run(['g++-16',str(OUT/'global_diagnostics.cpp'),'-Isrc/bench','-o',str(OUT/'global_diagnostics'),*base.FLAGS],check=True)
    with (OUT/'global_warning_details.log').open('w') as log:
        subprocess.run([str(OUT/'global_diagnostics'),str(OUT/'global_warning_inputs.dat'),str(OUT/'global_warning_results.dat')],stderr=log,check=True)
    for r in base.rows(OUT/'global_warning_results.dat'):
        assert r[1]==rr[int(r[0])][1] and r[3]==rr[int(r[0])][3]


def finalize():
    manifest=json.loads((OUT/'manifest.json').read_text())
    for name,digest in manifest['sha256'].items():
        if name.startswith(('src/','table/')):
            assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest, 'Numerical source changed during benchmark'
    manifest['final_script_sha256']={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in ['py/compare_all_optical_depth.py','py/compare_qagp.py']}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))


if __name__ == '__main__':
    actions = dict(setup=setup, run=run, analyze=analyze, refine=refine, audit=audit, audit_analyze=audit_analyze,
                   velocity=base.velocity, originals=base.originals, validate_ports=validate_ports,
                   standalone=standalone, diagnose_warnings=diagnose_warnings, finalize=finalize)
    for name in sys.argv[1:] or ['setup', 'run', 'velocity', 'analyze', 'refine', 'audit', 'audit_analyze', 'originals', 'validate_ports', 'standalone', 'diagnose_warnings', 'finalize']:
        actions[name]()
