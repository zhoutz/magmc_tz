import numpy as np
from scipy.integrate import solve_bvp

# documentation:
# https://typst.app/project/rNtvdGP93moPZ5T85QYooE


def solve_twisted_dipole(
    delta_phi: float,
    mu_out: np.ndarray,
    eps=1e-4,
    n_mesh=int(1e3),
    n_cont=100,
    tol=1e-7,
    max_nodes=int(1e5),
    verbose=False,
):
    """
    Solve the self-similar twisted dipole angular BVP for a target twist angle.

    Input
    -----
    delta_phi : float
        Target Delta phi_NS in radians. This routine assumes delta_phi >= 0.
    mu_out : np.ndarray
        Grid points where f and f' are returned. mu = cos(theta), 0 <= mu <= 1.
    eps : float
        Solve on [0, 1 - eps] and use the polar expansion near mu = 1.
        Keep this away from machine precision: too small a cutoff causes
        roundoff and excessive mesh refinement near the pole.
    n_mesh : int
        Initial collocation mesh size.
    n_cont : int
        Number of continuation steps in delta_phi.
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

    if delta_phi < 0:
        raise ValueError("This routine assumes delta_phi >= 0.")

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
    x = np.linspace(0.0, mu_max, n_mesh)

    # Initial guess: dipole field.
    f0 = (1.0 - x) * (1.0 + x)
    g0 = -2.0 * x
    J0 = x.copy()

    y_init = np.vstack([f0, g0, J0])
    n_cont = max(
        n_cont, 4, int(np.ceil(delta_phi / 0.05))
    )  # at least 4 steps, or more for large delta_phi
    # Match the first continuation target, not the final twist angle.
    pA_init = np.array([1.0, delta_phi / (2.0 * n_cont)])

    for istep, target_delta_phi in enumerate(
        np.linspace(delta_phi / n_cont, delta_phi, n_cont)
    ):

        def fun(mu, y, pars):
            p, A = pars

            # Avoid invalid fractional powers from tiny negative numerical noise.
            f_pos = np.maximum(y[0], 1e-300)
            denom = (1.0 - mu) * (1.0 + mu)

            df = y[1]
            dg = -p * (p + 1.0) * (y[0] + A * A * f_pos ** (1.0 + 2.0 / p)) / denom
            dJ = f_pos ** (1.0 / p) / denom

            return np.vstack([df, dg, dJ])

        def bc(ya, yb, pars, target_delta_phi=target_delta_phi):
            p, A = pars
            s = eps

            # Polar expansion at mu = 1 - eps.
            f_pole = 2.0 * s - 0.5 * p * (p + 1.0) * s**2
            fp_pole = -2.0 + p * (p + 1.0) * s

            # Integral tail from mu = 1 - eps to mu = 1.
            J_tail = 2.0 ** (1.0 / p - 1.0) * p * s ** (1.0 / p)

            return np.array(
                [
                    ya[1],  # f'(0) = 0
                    ya[2],  # J(0) = 0
                    yb[0] - f_pole,  # f(1 - eps)
                    yb[1] - fp_pole,  # f'(1 - eps)
                    2.0 * A * (yb[2] + J_tail) - target_delta_phi,
                ]
            )

        sol = solve_bvp(
            fun,
            bc,
            x,
            y_init,
            p=pA_init,
            tol=tol,
            max_nodes=max_nodes,
            verbose=0,
        )

        if verbose:
            print(
                f"step {istep + 1}/{n_cont}, "
                f"target={target_delta_phi:.6g}, success={sol.success}, "
                f"p={sol.p[0]:.12g}, A={sol.p[1]:.12g}, "
                f"nodes={sol.x.size}, max_residual={np.max(sol.rms_residuals):.3g}"
            )

        if not sol.success:
            raise RuntimeError(
                f"Continuation step {istep + 1}/{n_cont}, "
                f"target={target_delta_phi:.6g}: {sol.message} "
                f"nodes={sol.x.size}, "
                f"max_residual={np.max(sol.rms_residuals):.3g}, eps={eps:.3g}. "
            )

        # Retain the adaptive mesh as well as the converged solution.
        x = sol.x
        y_init = sol.y
        pA_init = sol.p

    p, A = sol.p
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
    target_delta_phi = 3.0  # radians
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
