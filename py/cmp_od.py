import sys

import numpy as np

name = sys.argv[1]

data = np.loadtxt(f"output/{name}.txt")[:, -1]
truth = np.loadtxt("table/bench_od.txt")[:, -1]

rel_error = np.abs(data - truth) / np.abs(truth)
max_rel_error = np.max(rel_error)
mean_rel_error = np.mean(rel_error)

print(f"----- dataset '{name}' -----")
print(f"Max relative error: {max_rel_error:.2e}")
print(f"Mean relative error: {mean_rel_error:.2e}")
