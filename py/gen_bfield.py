import numpy as np
from scipy.integrate import solve_bvp

# documentation:
# https://typst.app/project/rNtvdGP93moPZ5T85QYooE


def solve_twisted_dipole(
    delta_phi: float,
    mu_out: np.ndarray,
    eps=1e-5,
    n_mesh=int(1e3),
    n_cont=20,
    tol=1e-8,
    max_nodes=int(1e5),
    verbose=False,
):
    """
    Solve the self-similar twisted dipole angular BVP for a target twist angle.

    Input
    -----
    delta_phi : float
        Target Delta phi_NS in radians, 0 <= delta_phi < pi. The pi endpoint
        is the singular split-monopole limit p -> 0.
    mu_out : np.ndarray
        Grid points where f and f' are returned. mu = cos(theta), 0 <= mu <= 1.
    eps : float
        Cutoff for the domain [0, 1 - eps]. Solve on [0, 1 - eps] and use the polar expansion near mu = 1.
        Keep this away from machine precision: too small a cutoff causes
        roundoff and excessive mesh refinement near the pole.
    n_mesh : int
        Initial collocation mesh size.
    n_cont : int
        Minimum number of continuation steps. Steps shrink near delta_phi=pi.
    tol : float
        solve_bvp tolerance.
    max_nodes : int
        Maximum number of nodes for solve_bvp.
    verbose : bool
        Print progress information.

    Return
    ------
    dict with keys:
        p, A, C, mu, f, fp
    """

    if not np.isfinite(delta_phi) or not 0.0 <= delta_phi < np.pi:
        raise ValueError(
            "delta_phi must satisfy 0 <= delta_phi < pi on the twisted-dipole "
            "branch; delta_phi=pi is the singular p=0 split-monopole limit."
        )

    if not 0.0 < eps < 1.0 or 1.0 - eps == 1.0:
        raise ValueError("eps must satisfy 0 < eps < 1 and 1 - eps < 1.")

    mu_out = np.asarray(mu_out, dtype=np.float64)

    if np.any(mu_out < 0.0) or np.any(mu_out > 1.0):
        raise ValueError("mu_out must satisfy 0 <= mu <= 1.")

    # Exact dipole solution.
    if abs(delta_phi) < 1e-15:
        f = 1.0 - mu_out**2
        fp = -2.0 * mu_out
        return {
            "p": 1.0,
            "A": 0.0,
            "C": 0.0,
            "mu": mu_out,
            "f": f,
            "fp": fp,
        }

    mu_max = 1.0 - eps
    s_pole = 1.0 - mu_max
    x = np.linspace(0.0, mu_max, n_mesh)

    # Use uniform steps in -log(1 - delta_phi/pi): near the p -> 0 limit,
    # a fixed twist increment causes very large relative changes in p and A.
    stretch = -np.log1p(-delta_phi / np.pi)
    n_cont = max(n_cont, 4, int(np.ceil(np.pi * stretch / 0.05)))
    targets = -np.pi * np.expm1(np.linspace(0.0, -stretch, n_cont + 1)[1:])
    targets[-1] = delta_phi

    # Initial guess: dipole field.
    f0 = (1.0 - x) * (1.0 + x)
    g0 = -2.0 * x
    # Evolve Phi=2*A*J instead of J, which diverges as A -> 0 at large twist.
    phi0 = targets[0] * x
    y_init = np.vstack([f0, g0, phi0])
    # Solve for log(A) so tiny positive amplitudes remain well scaled.
    pars_init = np.array([1.0, np.log(targets[0] / 2.0)])

    for istep, target_delta_phi in enumerate(targets):

        def fun(mu, y, pars):
            p, log_A = pars

            # Avoid invalid fractional powers from tiny negative numerical noise.
            f_pos = np.maximum(y[0], 1e-300)
            denom = (1.0 - mu) * (1.0 + mu)
            # Evaluate A*f**(1/p) together, avoiding a huge power times tiny A.
            toroidal_ratio = np.exp(log_A + np.log(f_pos) / p)

            df = y[1]
            dg = -p * (p + 1.0) * (y[0] + f_pos * toroidal_ratio**2) / denom
            dphi = 2.0 * toroidal_ratio / denom

            return np.vstack([df, dg, dphi])

        def bc(ya, yb, pars, target_delta_phi=target_delta_phi):
            p, log_A = pars
            s = s_pole

            # Polar expansion at mu = 1 - eps.
            f_pole = 2.0 * s - 0.5 * p * (p + 1.0) * s**2
            fp_pole = -2.0 + p * (p + 1.0) * s

            # Phi tail = 2*A*J_tail = A*p*(2*s)**(1/p).
            phi_tail = p * np.exp(log_A + np.log(2.0 * s) / p)

            return np.array(
                [
                    ya[1],  # f'(0) = 0
                    ya[2],  # Phi(0) = 0
                    yb[0] - f_pole,  # f(1 - eps)
                    yb[1] - fp_pole,  # f'(1 - eps)
                    yb[2] + phi_tail - target_delta_phi,
                ]
            )

        sol = solve_bvp(
            fun,
            bc,
            x,
            y_init,
            p=pars_init,
            tol=tol,
            max_nodes=max_nodes,
            verbose=0,
        )

        if verbose:
            print(
                f"step {istep + 1}/{n_cont}, "
                f"target={target_delta_phi:.6g}, success={sol.success}, "
                f"p={sol.p[0]:.12g}, A={np.exp(sol.p[1]):.12g}, "
                f"nodes={sol.x.size}, max_residual={np.max(sol.rms_residuals):.3g}"
            )

        if not sol.success:
            raise RuntimeError(
                f"Continuation step {istep + 1}/{n_cont}, "
                f"target={target_delta_phi:.6g}: {sol.message} "
                f"nodes={sol.x.size}, "
                f"max_residual={np.max(sol.rms_residuals):.3g}, eps={eps:.3g}. "
            )

        # Keep the solution, but interpolate back to the base mesh each step.
        # solve_bvp only adds nodes; retaining every intermediate Newton mesh
        # accumulates unnecessary nodes and can exhaust max_nodes at large twist.
        y_init = sol.sol(x)
        pars_init = sol.p

    p, log_A = sol.p
    A = np.exp(log_A)
    C = A * A * p * (p + 1.0)

    # Evaluate on requested grid.
    mu_eval = np.minimum(mu_out, mu_max)
    y_eval = sol.sol(mu_eval)

    f = y_eval[0].copy()
    fp = y_eval[1].copy()

    # Fill the exact endpoint mu=1, or points extremely close to it, by expansion.
    mask = mu_out > mu_max
    if np.any(mask):
        s = 1.0 - mu_out[mask]
        f[mask] = 2.0 * s - 0.5 * p * (p + 1.0) * s**2
        fp[mask] = -2.0 + p * (p + 1.0) * s

    return {
        "p": float(p),
        "A": float(A),
        "C": float(C),
        "mu": mu_out,
        "f": f,
        "fp": fp,
    }


if __name__ == "__main__":
    target_delta_phi = 1.0  # radians
    mu_min = 0.0
    mu_max = 1.0
    mu_num = 10001

    mu_grid = np.linspace(mu_min, mu_max, mu_num)
    output_path = f"table/bfield_t{int(target_delta_phi * 10):02d}.txt"

    out = solve_twisted_dipole(
        delta_phi=target_delta_phi,
        mu_out=mu_grid,
        verbose=True,
    )

    print("p =", out["p"])
    print("A =", out["A"])
    print("C =", out["C"])
    print("f(0) =", out["f"][0])
    print("f'(0) =", out["fp"][0])
    print("f(1) =", out["f"][-1])
    print("f'(1) =", out["fp"][-1])

    with open(output_path, "w") as f:
        n = len(out["mu"])
        f.write(
            f"{target_delta_phi} {out['p']} {out['A']} {out['C']} {mu_min} {mu_max} {mu_num}\n"
        )
        f.writelines(f"{out['f'][i]:.16e} {out['fp'][i]:.16e}\n" for i in range(n))
    print(f"Table saved to {output_path}")

    if True:
        import matplotlib.pyplot as plt

        plt.plot(out["mu"], out["f"], label=r"f($\mu$)")
        plt.plot(out["mu"], out["fp"], label=r"f'($\mu$)")
        plt.xlabel(r"$\mu$")
        plt.ylabel("f, f'")
        plt.title(f"Twisted Dipole: $\\Delta\\phi$={target_delta_phi:.3g} rad")
        plt.legend()
        plt.grid()
        plt.savefig(f"table/bfield_t{int(target_delta_phi * 10):02d}.png", dpi=400)
        print(f"Plot saved to table/bfield_t{int(target_delta_phi * 10):02d}.png")
