// nonradial_ref.cpp — reference output for non-radial photon paths
//
// Uses the same algorithm as ref.cpp (3-component orbit + QAGS per step,
// max_step = 0.01r) but sweeps alpha0 in addition to the usual parameters.
// Output: output/nonradial_ref.txt
// Format: b0 muz alpha0 oi pol tau

#include "../bfield.hpp"
#include "../constants.hpp"
#include "../distribution.hpp"
#include "../dopr5.hpp"
#include "../photon.hpp"
#include "../qags.hpp"
#include "../solve_quadratic.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <print>

constexpr double M_star = 1.4;
constexpr double R_star = 10;
constexpr double rs     = M_star * schwarzschild_radius_of_sun_in_km;
constexpr double B_pole = 1e14;

using YVector = std::array<double, 3>;

struct PhotonEvolution {
    BField const    &bfield;
    Boltzmann const &fb;
    double3  n, e1, e2;
    double   omega_inf;
    Polarization pol;

    void operator()(double /*x*/, YVector const &y, YVector &dydx) const {
        auto [r, psi, alpha] = y;
        double f = std::sqrt(1 - rs / r);
        dydx[0] = f * std::cos(alpha);
        dydx[1] = std::sin(alpha) / r;
        dydx[2] = -std::sin(alpha) / (r * f) * (1 - 3 * rs / (2 * r));
    }

    double calc_dtaudl(double3 n_, double3 e1_, double3 e2_,
                       double r, double psi, double alpha,
                       double omega_inf_, Polarization pol_) const {
        double3 r_hat = std::cos(psi) * e1_ + std::sin(psi) * e2_;
        double  muz_  = r_hat.z;
        double3 B_vec = bfield.calc_B(r, muz_);
        double  B     = B_vec.length();
        double3 b     = B_vec / B;
        double  omega_c = B_to_omega * B;
        double  omega   = omega_inf_ / std::sqrt(1 - rs / r);
        double  x       = omega_c / omega;
        double  rho     = std::hypot(r_hat.x, r_hat.y);
        double3 theta_hat = rho > 0
            ? double3{r_hat.x*muz_/rho, r_hat.y*muz_/rho, -rho}
            : double3{std::copysign(1.0, muz_), 0, 0};
        double3 phi_hat = rho > 0
            ? double3{-r_hat.y/rho, r_hat.x/rho, 0}
            : double3{0, 1, 0};
        double mu_in = b.x*std::cos(alpha)
            + std::sin(alpha)*(b.y*dot(n_, phi_hat) - b.z*dot(n_, theta_hat));
        std::array<double, 2> betas;
        if (!solve_quadratic(x*x + mu_in*mu_in, -2*mu_in, (1+x)*(1-x), betas)) return 0;
        double ret = 0;
        for (double beta : betas) {
            double fv = fb.f(beta);
            if (fv == 0) continue;
            double mup_in = (mu_in - beta) / (1 - beta*mu_in);
            double esq    = (pol_ == Polarization::E) ? 0.5 : 0.5*mup_in*mup_in;
            ret += fv * esq * (1 - beta*mu_in)*(1 - beta*mu_in)
                   * (1 + beta)*(1 - beta) / std::abs(mu_in - beta);
        }
        ret *= (bfield.p + 1) * pi * bfield.Bphi_over_Btheta(muz_)
               / (std::abs(fb.b_bar()) * r);
        if (!std::isfinite(ret))
            throw std::runtime_error("Non-finite dtaudl encountered");
        return ret;
    }
};

double event_escape(double /*x*/, YVector const &y) { return y[0] - 10000; }

BField bfield("table/bfield_t10.txt", B_pole, R_star);

double total_optical_depth(double b0, double muz_init, double alpha0,
                           double oi, Polarization pol) {
    Boltzmann fb(b0);
    double3 r_hat{std::sqrt(1 - muz_init*muz_init), 0, muz_init};
    double3 n{0, 1, 0};

    PhotonEvolution pe{
        .bfield    = bfield,
        .fb        = fb,
        .n         = n,
        .e1        = r_hat,
        .e2        = cross(n, r_hat),
        .omega_inf = oi,
        .pol       = pol,
    };

    StepperDopr5<3, PhotonEvolution> stepper(pe, 1e-6, 1e-6);
    stepper.add_event(&event_escape);

    stepper.add_event([&](double, YVector const &y) {
        auto [r, psi, alpha] = y;
        double3 rh  = std::cos(psi)*pe.e1 + std::sin(psi)*pe.e2;
        double  muz_ = rh.z;
        double3 Bv  = bfield.calc_B(r, muz_);
        double  B   = Bv.length();
        double3 b   = Bv / B;
        double  xv  = B_to_omega*B / (oi/std::sqrt(1 - rs/r));
        double  rho = std::hypot(rh.x, rh.y);
        double3 th  = rho > 0
            ? double3{rh.x*muz_/rho, rh.y*muz_/rho, -rho}
            : double3{std::copysign(1.0, muz_), 0, 0};
        double3 ph  = rho > 0
            ? double3{-rh.y/rho, rh.x/rho, 0}
            : double3{0, 1, 0};
        double mu_in = b.x*std::cos(alpha)
            + std::sin(alpha)*(b.y*dot(n, ph) - b.z*dot(n, th));
        return std::fma(xv, xv, (mu_in - 1)*(mu_in + 1));
    });

    stepper.add_event([&](double, YVector const &y) {
        double r = y[0];
        return B_to_omega * bfield.calc_B(r, muz_init).length()
               * std::sqrt(1 - rs/r) / oi - 1;
    });

    const bool stop_at_beta_zero = b0 * bfield.calc_B(R_star, muz_init).x <= 0;
    stepper.init(0.0, 1e-3*R_star, {R_star, 0.0, alpha0});
    double tau = 0;

    while (true) {
        stepper.do_step(1e-2 * stepper.y_old[0]);
        stepper.prepare_dense();
        int event_id = stepper.detect_event();

        double dtau = qags(
            [&](double s) {
                auto [r, psi, alpha] = stepper.dense_out(s);
                return pe.calc_dtaudl(pe.n, pe.e1, pe.e2, r, psi, alpha,
                                      pe.omega_inf, pe.pol);
            },
            stepper.x_old, stepper.x_old + stepper.h_old, 1e-6, 1e-6);

        tau += dtau;

        if (event_id == 0 || event_id == 1
            || (event_id == 2 && stop_at_beta_zero)) break;
        stepper.update_old();
    }
    return tau;
}

int main() {
    std::FILE *fp = std::fopen("output/nonradial_ref.txt", "w");
    if (!fp) { std::perror("output/nonradial_ref.txt"); return 1; }

    for (double b0   : {-0.3, -0.6})
    for (double muz  : {0.0, 0.3, 0.6, 0.9})
    for (double a0   : {0.0, pi/6, pi/3})
    for (double oi   : {0.1, 1., 10.})
    for (Polarization pol : {Polarization::E, Polarization::O}) {
        double tau = total_optical_depth(b0, muz, a0, oi, pol);
        std::println(fp, "{:.2f} {:.2f} {:.6f} {:.2f} {} {:.16e}",
                     b0, muz, a0, oi,
                     (pol == Polarization::E ? 1 : 0), tau);
    }

    std::fclose(fp);
    return 0;
}
