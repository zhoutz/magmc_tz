#pragma once

// Outward radial rays only. mu and magnetic colatitude stay constant.
#include "photon_evolution.hpp"
#include "roots.hpp"
#include <algorithm>
#include <limits>

struct RadialResonance {
  Boltzmann const &fb;
  double mu, q, surface_x, prefactor;
  Polarization pol;

  RadialResonance(BField const &field, Boltzmann const &distribution,
                  double muz, double energy, Polarization mode)
      : fb(distribution), pol(mode) {
    auto B = field.calc_B(R_star, muz);
    mu = B.x / B.length();
    q = field.p + 2;
    surface_x = B.length() * B_to_omega / energy;
    prefactor = (field.p + 1) * pi * field.Bphi_over_Btheta(muz) / std::abs(fb.b_bar());
  }
  double lapse(double r) const { return std::sqrt(1 - rs / r); }
  double Q(double r) const { return q - 0.5 * rs / (r - rs); }
  double x(double r) const { return surface_x * std::pow(R_star / r, q) * lapse(r); }
  double g(double beta) const { return (1 - beta * mu) / std::sqrt((1-beta)*(1+beta)); }
  double beta(double r) const {
    double X = x(r);
    if (X <= 1) return 0;
    // Rationalized negative branch, stable near beta=0.
    double c = (X-1)*(X+1);
    return -c / (mu + X * std::sqrt(c + mu*mu));
  }
  double radius_for_x(double X) const {
    return zriddr([&](double r) { return std::log(x(r)/X); }, R_star, 10000., 2e-12);
  }
  double dbeta_dl(double r, double b) const {
    return -Q(r)*lapse(r)/r * (1-b*mu)*(1-b*b)/(b-mu);
  }
  double dtau_dr(double r) const {
    const double b = beta(r), f = fb.f(b);
    if (f == 0) return 0;
    const double m = (mu-b)/(1-b*mu);
    const double P = pol == Polarization::E ? .5 : .5*m*m;
    return prefactor/r/lapse(r) * f*P*(1-b*mu)*(1-b*mu)*(1-b*b)/std::abs(mu-b);
  }
};

// FT07 sections 3.3--3.4, equations (31), (38), (39).
// The paper leaves distribution edges and the beta=0 constant empirical.
// Here edges enclose the central 99.8%, and the centre is the median.
// GR adaptation: d ln(omega_c/omega)/dl = -Q L/r, not just d ln B/dl.
struct FT07StepControl {
  RadialResonance const &ray;
  double low, centre, high;
  double resolution = .01;

  double operator()(double r) const {
    double cap = r/10;
    double b = ray.beta(r);
    if (b == 0) return cap;
    double speed = std::abs(ray.dbeta_dl(r, b));
    if (b >= low && b <= high) {
      double db = resolution*(1-b*b)*std::abs(b); // Eq. 39
      if (db == 0) db = resolution*std::min(centre-low, high-centre);
      db = std::min(db, resolution*std::abs(ray.mu-b)); // Eq. 38
      cap = std::min(cap, db/speed);                 // Eq. 31
    } else if (b < low) {
      // Resolve an approaching edge even if BOTH ends of the large step
      // would lie outside the resonant interval. Monotonic radial x makes
      // this a bracketed geometric test rather than an opacity test.
      double probe = r + cap; // dr/dl <= 1: conservative look-ahead
      if (ray.beta(probe) >= low) {
        double db = std::min(low-b + resolution*(centre-low), resolution*std::abs(ray.mu-b));
        cap = std::min(cap, db/speed);
      }
    }
    return cap;
  }
};

// A cheap full-distribution alternative: cap the change in log(x) according
// to the cold equatorial layer's O(beta0^2) width. The edge is used only to
// relax the cap in the far tail; it never truncates the optical depth.
struct ThermalStepControl {
  RadialResonance const &ray;
  double tail_beta;
  double resolution=.01;
  double operator()(double r) const {
    double X=ray.x(r), cap=r/10;
    if(X>1) {
      double width=resolution*ray.fb.b0*ray.fb.b0;
      cap=std::min(cap,(width+std::max(0.,std::log(X)-std::log(ray.g(tail_beta))))*
                   r/(ray.Q(r)*ray.lapse(r)));
    }
    return cap;
  }
};
