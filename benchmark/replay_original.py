#!/usr/bin/env python3
"""Reproduce the original implementation from the pre-fix Git commit."""
from pathlib import Path
import os
import subprocess

commit = '2789b0e1cc075fd19241d75e5adff9771880c54e'
root = Path(__file__).resolve().parents[1]
folder = root / 'build' / 'original_snapshot'
folder.mkdir(parents=True, exist_ok=True)
files = subprocess.check_output(['git', 'ls-tree', '--name-only', commit + ':ode'], cwd=root, text=True).splitlines()
for name in files:
    if name.endswith(('.cpp', '.hpp')):
        (folder / name).write_bytes(subprocess.check_output(['git', 'show', commit + ':ode/' + name], cwd=root))
(folder / 'probe.cpp').write_text(r'''#define main original_main
#include "main.cpp"
#undef main
#include <iostream>
#include <iomanip>
int main() {
  std::cout << std::setprecision(17) << "beta0,muz,pol,tau\n";
  for (double beta0 : {-.2, -.005}) for (double muz : {-.2, 0.})
    for (auto mode : {Polarization::O, Polarization::E}) {
      Boltzmann distribution(beta0);
      Ran random(1);
      double st = std::sqrt(1-muz*muz);
      PhotonEvolution physics{bfield, distribution, random, {0,1,0}, {st,0,muz}, {muz,0,-st},
                              {10,0,0,-1}, 1, mode};
      StepperDopr5<4,PhotonEvolution> stepper(physics,1e-6,1e-6);
      stepper.init(0,.01,physics.r_psi_alpha_tau);
      for (int k=0; k<1000000 && stepper.y_old[0]<10000; ++k) {
        stepper.do_step(); stepper.update_old();
      }
      if (stepper.y_old[0]<10000) throw std::runtime_error("step limit");
      std::cout << beta0 << ',' << muz << ',' << (mode==Polarization::O?'O':'E') << ','
                << stepper.y_old[3]+1 << '\n';
    }
}
''')
subprocess.run([os.environ.get('CXX', 'g++-16'), '-std=c++23', '-O3', str(folder/'probe.cpp'), '-o', str(folder/'probe')], cwd=root, check=True)
result = subprocess.check_output([str(folder/'probe')], cwd=root, text=True)
(root/'benchmark'/'original.csv').write_text(result)
print(result, end='')
