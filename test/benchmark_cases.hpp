#pragma once
#include "../ode/transport_physics.hpp"
#include <string>
#include <vector>

struct Case {
  std::string name;
  double b0;
  Photon photon;
  bool radial;
};

Photon make_photon(double muz, double alpha, double azimuth, double energy, Polarization mode,
                   double radius = R_star) {
  double st = std::sqrt((1 - muz) * (1 + muz));
  double3 e1{st, 0, muz}, theta{muz, 0, -st}, phi{0, 1, 0};
  double3 e2 = std::cos(azimuth) * theta + std::sin(azimuth) * phi;
  return {cross(e1, e2), e1, e2, radius, 0, alpha, energy, mode};
}

std::vector<Case> cases(bool quick) {
  std::vector<Case> result;
  auto add = [&](std::string label, double b0, double muz, double energy, double alpha = 0,
                 double azimuth = 0, double radius = R_star) {
    for (auto mode : {Polarization::O, Polarization::E}) {
      result.push_back({label + (mode == Polarization::O ? "_O" : "_E"), b0,
                        make_photon(muz, alpha, azimuth, energy, mode, radius),
                        std::abs(std::sin(alpha)) < 1e-14});
    }
  };
  add("review", -.2, -.2, 1);
  if (quick) {
    add("cold_equator", -.005, 0, 1);
    add("ultracold_south", -.001, -.2, 10);
    add("hot_north", -.75, .8, .1);
    add("oblique", -.2, -.2, 1, 1.2, 0.7);
    add("inward", -.05, -.2, 1, pi, 0, 300);
    return result;
  }
  int i = 0;
  for (double b0 : {-.75, -.2, -.05, -.005})
    for (double muz : {-.8, -.2, 0., .2, .8})
      for (double energy : {.1, 1., 10.}) add("grid" + std::to_string(i++), b0, muz, energy);
  for (double muz : {-.2, 0., .2}) add("ultracold" + std::to_string(i++), -.001, muz, 1);
  for (double muz : {-1., 1.}) add("axis" + std::to_string(i++), -.2, muz, 1);
  for (double b0 : {-.75, -.05}) {
    add("inward_south" + std::to_string(i++), b0, -.2, 1, pi, 0, 300);
    add("inward_north" + std::to_string(i++), b0, .2, 1, pi, 0, 300);
    for (double az : {0., 1.2}) {
      add("oblique_south" + std::to_string(i++), b0, -.2, 1, 1.2, az);
      add("oblique_north" + std::to_string(i++), b0, .2, 10, .9, az);
      add("turning" + std::to_string(i++), b0, -.2, 1, 2.4, az, 140);
    }
  }
  add("positive_carriers", .2, .2, 1);
  // Starts inside/near the layer; phase and initial proposed step differ.
  add("inside_layer", -.2, -.2, 1, 0, 0, 90);
  // Choose actual resonance positions, rather than approximate radii, to test
  // truncated layers and launch-phase dependence against the velocity reference.
  BField field("table/bfield_t10.txt", B_pole, R_star);
  for (double muz : {-.2, 0.}) {
    double3 B = field.calc_B(R_star, muz);
    double mu = B.x / B.length();
    for (double beta : {-.1, -.001}) {
      double target = (1 - beta * mu) / std::sqrt(1 - beta * beta);
      double lo = R_star, hi = 300;
      for (int k = 0; k < 64; ++k) {
        double r = .5 * (lo + hi);
        double x = B_to_omega * field.calc_B(r, muz).length() * std::sqrt(1-rs/r);
        if (x > target) lo = r; else hi = r;
      }
      add("partial_layer" + std::to_string(i++), beta == -.1 ? -.2 : -.001,
          muz, 1, 0, 0, .5 * (lo + hi));
    }
  }
  add("hard_cold", -.005, 0, 100);
  add("inward_surface", -.2, .2, 10, 2.99, .4, 30);
  return result;
}

