from itertools import product

import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar


class Bfield:
    def __init__(self, fname, B_pole, R_star):
        self.B_pole = B_pole
        self.R_star = R_star
        with open(fname, "r") as f:
            header = f.readline().split()
            self.Delta_phi = float(header[0])
            self.p = float(header[1])
            self.A = float(header[2])
            self.C = float(header[3])
            self.mu_min = float(header[4])
            self.mu_max = float(header[5])
            self.mu_num = int(header[6])
            self.f = np.zeros(self.mu_num)
            self.fp = np.zeros(self.mu_num)
            for i in range(self.mu_num):
                line = f.readline().split()
                self.f[i] = float(line[0])
                self.fp[i] = float(line[1])

    def _interpolate(self, mu):
        mu_abs = np.abs(mu)
        idx_float = (
            (mu_abs - self.mu_min) / (self.mu_max - self.mu_min) * (self.mu_num - 1)
        )
        i = int(np.floor(idx_float))
        i = np.clip(i, 0, self.mu_num - 2)
        a = idx_float - i
        fval = (1 - a) * self.f[i] + a * self.f[i + 1]
        fpval = (1 - a) * self.fp[i] + a * self.fp[i + 1]
        return fval, fpval

    def calc_B(self, r, mu):
        mu_sign = 1.0 if mu >= 0 else -1.0
        fval, fpval = self._interpolate(mu)
        fpval *= mu_sign
        sth = np.sqrt(1 - mu * mu)
        Br = -fpval
        Bth = self.p * fval / sth if sth != 0 else 0.0
        Bph = self.A * np.power(fval, 1.0 / self.p) * Bth
        scale = 0.5 * self.B_pole * np.power(self.R_star / r, 2 + self.p)
        return scale * Br, scale * Bth, scale * Bph

    def Bphi_over_Btheta(self, mu):
        fval, _ = self._interpolate(mu)
        return self.A * np.power(fval, 1.0 / self.p)


class Boltzmann:
    def __init__(self, b0):
        from scipy.special import k1e

        if not np.isfinite(b0) or np.abs(b0) >= 1 or b0 == 0:
            raise ValueError("Invalid b0 value for Boltzmann distribution")

        self.b0 = b0

        if b0 < 0:
            self.b_min = -1.0
            self.b_max = 0.0
        else:
            self.b_min = 0.0
            self.b_max = 1.0

        s0 = np.sqrt((1 - b0) * (1 + b0))
        self.a = s0 * (1 + s0) / (b0 * b0)
        self.exp_k1 = k1e(self.a)  # Scaled modified Bessel function K1
        self._b_bar = np.copysign(1.0 / (self.a * self.exp_k1), b0)

    def f(self, b):
        if not self.b_min < b < self.b_max:
            return 0.0

        s = np.sqrt((1 - b) * (1 + b))
        gm1 = b * b / (s * (1 + s))
        return np.exp(-self.a * gm1) / (self.exp_k1 * s * s * s)

    def b_bar(self):
        return self._b_bar


R_star = 10
B_pole = 1e14
rlo, rhi = R_star, 10000
bfield = Bfield("table/bfield_t10.txt", B_pole=B_pole, R_star=R_star)
q = bfield.p + 2
schwarzschild_radius_of_sun_in_km = 2.9532501
rs = 1.4 * schwarzschild_radius_of_sun_in_km
B_to_omega = 1.157676359638892e-11


def calc_optical_depth(b0, muz, omega_inf, pol):
    boltzmann = Boltzmann(b0=b0)
    Bvec = bfield.calc_B(R_star, muz)
    B = np.sqrt(Bvec[0] ** 2 + Bvec[1] ** 2 + Bvec[2] ** 2)
    mu = Bvec[0] / B
    omega_c = B * B_to_omega

    def integral(beta):
        g = (1 - beta * mu) / np.sqrt(1 - beta**2)

        def x(r):
            lapse = np.sqrt(1 - rs / r)
            return omega_c * lapse / omega_inf * (R_star / r) ** q

        if not x(rhi) < g < x(rlo):
            return 0.0

        r = root_scalar(
            lambda r: x(r) - g,
            bracket=[rlo, rhi],
            method="brentq",
            xtol=1e-10,
            rtol=1e-10,
        ).root
        L = np.sqrt(1 - rs / r)
        Q = q - 0.5 * rs / (r - rs)

        if pol == "E":
            P = 0.5
        elif pol == "O":
            P = 0.5 * ((mu - beta) / (1 - beta * mu)) ** 2
        else:
            raise ValueError("Invalid polarization type. Use 'E' or 'O'.")

        prefactor = (
            (bfield.p + 1)
            * np.pi
            * bfield.Bphi_over_Btheta(muz)
            / np.abs(boltzmann.b_bar())
        )
        body = boltzmann.f(beta) * P * (1 - beta * mu) / (L * Q)

        return prefactor * body

    return quad(integral, boltzmann.b_min, boltzmann.b_max, epsabs=0, epsrel=1e-11)[0]


output_path = "table/bench_od.txt"
with open(output_path, "w") as f:
    # f.write("# b0 muz pol omega_inf tau\n")
    for b0, muz, oi, pol in product(
        np.linspace(-0.1, -0.9, 9),
        np.linspace(0, 0.9, 10),
        [0.01, 0.1, 1, 10, 100],
        ["E", "O"],
    ):
        tau = calc_optical_depth(b0, muz, oi, pol)
        print(
            f"b0={b0:.2f}, muz={muz:.2f}, pol={pol}, omega_inf={oi:.2f}, tau={tau:.16e}"
        )
        f.write(f"{b0:.2f} {muz:.2f} {oi:.2f} {0 if pol == 'O' else 1} {tau:.16e}\n")
    print(f"Results written to {output_path}")
