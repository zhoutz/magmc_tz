#pragma once

#include "photon.hpp"
#include "ran.hpp"

inline Photon init07(Ran &ran, double R_star, Polarization pol) {
  double3 e1 = ran.point_on_unit_sphere();
  double3 n = ran.unit_perp_to(e1);
  double3 e2 = cross(n, e1);
  return Photon{.n = n,
                .e1 = e1,
                .e2 = e2,
                .r = R_star,
                .psi = 0.0,
                .alpha = 0.0,
                .omega_inf = 1.0,
                .pol = pol};
}
