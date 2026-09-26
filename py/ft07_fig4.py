#!/usr/bin/env python3
"""FT07 Fig. 4(b, e, f), beta0=0.75 and Delta_phi=1.

From the repository root:
  build/ft07-venv/bin/python py/ft07_fig4.py
  build/ft07-venv/bin/python py/ft07_fig4.py --check

Dependencies: numpy, matplotlib. No scipy is required.
FT07 sections 3.6 and Appendix A: number response -> Planck convolution ->
9-degree observer cone -> rotational average -> energy per log10 frequency.
The input blackbody has unit energy integral over log10(omega/omega_bb).
All estimates use the emitted count, retaining the loss to surface absorption.
The shaded bands are one Monte Carlo standard error (bin approximation), not
an estimate of discretization error or of an unsampled high-energy tail.
"""

import argparse
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "build" / "ft07-mpl-cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / "build" / "ft07-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.legendre import leggauss

ZETA3 = 1.2020569031595943
MEAN_X = np.pi**4 / (30 * ZETA3)


def blackbody_number(log_x):
    """Unit photon-number density per dex; stable in both Planck tails."""
    log_x = np.asarray(log_x)
    x = 10.0 ** np.clip(log_x, -100, np.log10(700.0))
    value = np.log(10) / (2 * ZETA3) * x**3 * np.exp(-x) / (-np.expm1(-x))
    return np.where(log_x > np.log10(700.0), 0.0, value)


def read_photons(path):
    metadata = {}
    with path.open() as stream:
        for line in stream:
            if line.startswith("#") and "=" in line:
                key, value = line[1:].split("=", 1)
                metadata[key.strip()] = value.strip()
    if metadata.get("complete") != "1":
        raise ValueError("Input is incomplete: missing '# complete = 1' footer")
    n = int(metadata["n_emitted"])
    escaped = int(metadata["n_escaped"])
    absorbed = int(metadata["n_absorbed"])
    if n <= 1 or escaped + absorbed != n:
        raise ValueError("Invalid emitted/escaped/absorbed counts")
    data = np.loadtxt(path, ndmin=2)
    if data.shape != (escaped, 5) or not np.all(np.isfinite(data)):
        raise ValueError(
            "Photon rows do not match metadata or contain nonfinite values"
        )
    if np.any(data[:, 0] <= 0) or np.any(np.abs(data[:, 1]) > 1):
        raise ValueError("Invalid photon energy or direction cosine")
    for key, expected in (("beta0", 0.75), ("delta_phi", 1.0), ("rs_km", 0.0)):
        if not np.isclose(float(metadata[key]), expected):
            raise ValueError(f"This figure requires {key}={expected}")
    return data, metadata


def make_response(data, metadata, ds, n_mu):
    s = np.log10(data[:, 0] / float(metadata["omega_in_keV"]))
    # Integer-multiple centers put the unscattered delta exactly at s=0.
    lo = int(np.floor(s.min() / ds))
    hi = int(np.ceil(s.max() / ds))
    centers = np.arange(lo, hi + 1) * ds
    s_edges = np.arange(lo, hi + 2) * ds - ds / 2
    mu_edges = np.linspace(-1, 1, n_mu + 1)
    counts = np.histogram2d(s, data[:, 1], bins=(s_edges, mu_edges))[0]
    if int(counts.sum()) != len(data):
        raise ValueError("Histogram dropped photons")
    response = counts / int(metadata["n_emitted"]) / ds / np.diff(mu_edges)[None, :]
    return counts, response, centers, s_edges, mu_edges


def cone_weights(mu_edges, mu_spin, mu_los, beam_deg=9.0, n_phase=512, order=64):
    """Angular-bin probabilities for a uniform cone averaged over spin phase.

    At magnetic cosine mu, the fraction of azimuth inside the cone is
    acos((cos(beam)-mu*v)/sqrt((1-mu^2)*(1-v^2)))/pi, clipped to [0,1],
    where v is the cone-axis magnetic cosine. Integrate this over each mu bin
    and divide by 1-cos(beam). This is FT07 eqs. (52)-(56) with the solid-angle
    normalization made explicit. No stochastic observer resampling is used.
    """
    c = np.cos(np.deg2rad(beam_deg))
    phase = 2 * np.pi * (np.arange(n_phase) + 0.5) / n_phase
    amplitude = np.sqrt((1 - mu_spin**2) * (1 - mu_los**2))
    observers = mu_spin * mu_los + amplitude * np.cos(phase)
    if amplitude == 0:
        observers = observers[:1]
    nodes, weights = leggauss(order)
    half = np.diff(mu_edges) / 2
    mu = (mu_edges[:-1] + half)[:, None] + half[:, None] * nodes
    integral = np.zeros(len(half))
    for v in observers:
        if abs(v) > 1 - 1e-14:
            # Exact polar cap overlap; avoids quadrature over a step function.
            left, right = (c, 1.0) if v > 0 else (-1.0, -c)
            integral += np.maximum(
                0, np.minimum(mu_edges[1:], right) - np.maximum(mu_edges[:-1], left)
            )
        else:
            z = (c - mu * v) / np.sqrt((1 - mu**2) * (1 - v**2))
            fraction = np.arccos(np.clip(z, -1, 1)) / np.pi
            integral += half * (fraction @ weights)
    result = integral / (len(observers) * (1 - c))
    normalization_error = float(result.sum() - 1)
    if abs(normalization_error) > 2e-3:
        raise ValueError(
            "Observer-cone quadrature is underresolved; increase --angle-order"
        )
    # Enforce exact isotropic normalization after the small quadrature error.
    return result / result.sum(), normalization_error


def spectrum_and_error(counts, kernel, angular_weights, mu_edges, n, x):
    # An individual photon in bin (i,j) contributes
    # kernel(t,s_i) * angular_weights[j]/dmu[j] to the angular number density.
    # Absorbed photons have zero contribution but are included in n.
    a = angular_weights / np.diff(mu_edges)
    first = kernel @ (counts @ a)
    second = kernel**2 @ (counts @ a**2)
    factor = 2 * x / MEAN_X
    spectrum = factor * first / n
    variance_of_mean = np.maximum(second - first**2 / n, 0) / (n * (n - 1))
    return spectrum, factor * np.sqrt(variance_of_mean)


def check_numerics():
    t = np.linspace(-6, 3, 20001)
    b = blackbody_number(t)
    np.testing.assert_allclose(np.trapezoid(b, t), 1, rtol=1e-8)
    np.testing.assert_allclose(np.trapezoid(10**t * b / MEAN_X, t), 1, rtol=1e-8)
    edges = np.linspace(-1, 1, 201)
    n = 200000
    counts = np.full((1, 200), n / 200)
    t = np.linspace(-1, 2, 151)
    x = 10**t
    # Independent analytic Planck energy shape: unscattered isotropic photons
    # must give this spectrum for every geometry, including polar observers.
    expected = 15 * np.log(10) / np.pi**4 * x**4 / np.expm1(x)
    for spin, los in ((0, 0), (1, 1), (1, -1), (1, 0.5), (0, 0.5)):
        w, _ = cone_weights(edges, spin, los)
        actual, _ = spectrum_and_error(
            counts, blackbody_number(t[:, None]), w, edges, n, x
        )
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-15)
    plus, _ = cone_weights(edges, 0, 0.5)
    minus, _ = cone_weights(edges, 0, -0.5)
    np.testing.assert_allclose(plus, minus, atol=1e-12)
    w1, _ = cone_weights(edges, 1, 0.5)
    w2, _ = cone_weights(edges, 1, -0.5)
    np.testing.assert_allclose(w1, w2[::-1], atol=1e-12)
    print(
        "Checks passed: Planck integrals, isotropic no-scattering limit, cone symmetries."
    )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=ROOT / "output/ft07_beta075_twist10_E.txt",
    )
    parser.add_argument("--output", type=Path, default=ROOT / "output/ft07_fig4.png")
    parser.add_argument("--ds", type=float, default=0.02)
    parser.add_argument("--n-mu", type=int, default=200)
    parser.add_argument("--beam-deg", type=float, default=9)
    parser.add_argument("--n-phase", type=int, default=512)
    parser.add_argument("--angle-order", type=int, default=64)
    parser.add_argument(
        "--check", action="store_true", help="Run analytic checks and exit"
    )
    args = parser.parse_args()
    if args.check:
        check_numerics()
        return
    if not (
        args.ds > 0
        and args.n_mu >= 2
        and 0 < args.beam_deg < 90
        and args.n_phase >= 4
        and args.angle_order >= 4
    ):
        parser.error("Invalid binning or angular quadrature parameters")
    data, metadata = read_photons(args.input)
    n = int(metadata["n_emitted"])
    counts, response, s, s_edges, mu_edges = make_response(
        data, metadata, args.ds, args.n_mu
    )
    t = np.arange(-1, 4.0001, 0.01)
    x = 10**t
    kernel = blackbody_number(t[:, None] - s[None, :])
    H = (kernel * np.diff(s_edges)[None, :]) @ response
    bb = x * blackbody_number(t) / MEAN_X
    curves, errors, observer_weights = [], [], []
    geometries = [(0.0, 0.0)]
    los_values = [1.0, 1 / np.sqrt(2), 0.0, -1 / np.sqrt(2), -1.0]
    geometries += [(spin, los) for spin in (1.0, 0.0) for los in los_values]
    quadrature_errors = []
    for spin, los in geometries:
        w, err = cone_weights(
            mu_edges, spin, los, args.beam_deg, args.n_phase, args.angle_order
        )
        y, se = spectrum_and_error(counts, kernel, w, mu_edges, n, x)
        np.testing.assert_allclose(y, 2 * x / MEAN_X * (H @ w), rtol=1e-12, atol=1e-15)
        curves.append(y)
        errors.append(se)
        observer_weights.append(w)
        quadrature_errors.append(err)
    curves, errors = np.array(curves), np.array(errors)
    fraction = float(
        np.sum(response * np.diff(s_edges)[:, None] * np.diff(mu_edges)[None, :])
    )
    np.testing.assert_allclose(fraction, len(data) / n, rtol=1e-12)

    plt.rcParams.update(
        {"font.size": 11, "axes.spines.top": True, "axes.spines.right": True}
    )
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharey=True)
    colors = ["#1b9e77", "#377eb8", "#222222", "#984ea3", "#d95f02"]
    styles = [(0, (7, 2, 7, 2, 1, 2)), "-.", "-", (0, (1, 2, 1, 2, 6, 2)), "--"]

    def draw(ax, index, color, style, label):
        y, se = curves[index], errors[index]
        ax.semilogy(
            t, np.where(y > 0, y, np.nan), color=color, ls=style, lw=1.5, label=label
        )
        ax.fill_between(
            t, np.maximum(y - se, 1e-12), y + se, color=color, alpha=0.12, linewidth=0
        )

    draw(axes[0], 0, colors[2], "-", r"$\mu_\Omega=\mu_{\rm los}=0$")
    labels = [r"$1$", r"$1/\sqrt{2}$", r"$0$", r"$-1/\sqrt{2}$", r"$-1$"]
    for i in range(5):
        draw(axes[1], i + 1, colors[i], styles[i], labels[i])
    # +/- sight lines coincide for an orthogonal rotator; show each once.
    for i, label in enumerate([r"$\pm1$", r"$\pm1/\sqrt{2}$", r"$0$"]):
        draw(axes[2], i + 6, colors[i], styles[i], label)
    titles = [
        r"(b) One twist: $\Delta\phi_{\rm N-S}=1$",
        r"(e) Aligned: $\mu_\Omega=1$",
        r"(f) Orthogonal: $\mu_\Omega=0$",
    ]
    for i, (ax, title) in enumerate(zip(axes, titles)):
        ax.semilogy(t, bb, ":", color="0.4", lw=1.7, label="Input blackbody")
        ax.set(
            xlim=(-1, 4),
            ylim=(0.02, 10),
            xlabel=r"$\log_{10}(\omega/\omega_{\rm bb})$",
            title=title,
        )
        ax.grid(alpha=0.18, which="major")
        ax.legend(
            fontsize=9, title=None if i == 0 else r"$\mu_{\rm los}$", frameon=False
        )
    axes[0].set_ylabel(r"$\omega F_\omega$ (normalized)")
    fig.suptitle(
        rf"FT07: $\beta_0=0.75$, $\Delta\phi_{{\rm N-S}}=1$; "
        rf"$N_{{\rm emitted}}={n:,}$, $\theta_{{\rm beam}}={args.beam_deg:g}^\circ$"
    )
    fig.text(
        0.5,
        0.015,
        "Shading: 1 Monte Carlo standard error. Input blackbody has unit area in log10 frequency.",
        ha="center",
        fontsize=9,
        color="0.35",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200)
    plt.close(fig)
    archive = args.output.with_suffix(".npz")
    np.savez_compressed(
        archive,
        log10_x=t,
        blackbody=bb,
        spectra=curves,
        standard_errors=errors,
        geometries=geometries,
        response=response,
        H=H,
        counts=counts,
        s_edges=s_edges,
        mu_edges=mu_edges,
        observer_weights=observer_weights,
        metadata_json=json.dumps(metadata),
        beam_deg=args.beam_deg,
        n_phase=args.n_phase,
        angle_order=args.angle_order,
    )
    summary = {
        "n_emitted": n,
        "n_escaped": len(data),
        "n_absorbed": int(metadata["n_absorbed"]),
        "escape_fraction": fraction,
        "response_integral": fraction,
        "max_cone_normalization_error": float(np.max(np.abs(quadrature_errors))),
        "photon_file": str(args.input),
        "figure": str(args.output),
        "archive": str(archive),
        "ds": args.ds,
        "dmu": 2 / args.n_mu,
        "beam_deg": args.beam_deg,
        "n_phase": args.n_phase,
        "angle_order": args.angle_order,
    }
    args.output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
