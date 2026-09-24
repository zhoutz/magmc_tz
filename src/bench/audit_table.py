"""Reproduce the radial-table discrepancy without executing its file-writing code."""
from pathlib import Path

source = Path("py/bench_od.py").read_text().split("output_path =")[0]
original = {}
# Only the definitions before the original script's output block are executed.
exec(compile(source, "py/bench_od.py", "exec"), original)
old = "return quad(integral, boltzmann.b_min, boltzmann.b_max, epsabs=0, epsrel=1e-11)[0]"
new = (
    "G = omega_c * np.sqrt(1-rs/R_star) / omega_inf\n"
    "    lower = -(G*G-1)/(mu+G*np.sqrt(G*G+mu*mu-1))\n"
    "    return quad(integral, lower, boltzmann.b_max, epsabs=0, epsrel=1e-11)[0]"
)
if old not in source:
    raise RuntimeError("The reference Python implementation changed; review the audit adapter")
corrected = {}
exec(compile(source.replace(old, new), "explicit_velocity_domain", "exec"), corrected)
with open("output/table_audit.txt", "w") as f:
    f.write("# b0 muz energy pol original_python explicit_domain_python\n")
    for b0 in [-.6, -.7, -.8, -.9]:
        for pol in ["E", "O"]:
            a = original["calc_optical_depth"](b0, .4, 100, pol)
            b = corrected["calc_optical_depth"](b0, .4, 100, pol)
            f.write(f"{b0} 0.4 100 {pol} {a:.17g} {b:.17g}\n")
print("Wrote output/table_audit.txt (original source/table unchanged)")
