import numpy as np
from scipy.integrate import quad

y_min = -4.0
y_max = 1.5
num_points = 551
output_path = "table/nbb_table.txt"


def integrand(x):
    if x == 0:
        return 0.0
    else:
        return x**2 / (np.exp(x) - 1)


y_vals = np.linspace(y_min, y_max, num_points)
cdf_vals = np.zeros_like(y_vals)

for i, y in enumerate(y_vals):
    x_val = 10**y
    integral_val, _ = quad(integrand, 0, x_val)
    cdf_vals[i] = integral_val

cdf_vals -= cdf_vals[0]
cdf_vals /= cdf_vals[-1]

saved_data = np.column_stack((y_vals, cdf_vals))
np.savetxt(output_path, saved_data, header="log10(E/kT) CDF")
print(f"Table saved to {output_path}")
