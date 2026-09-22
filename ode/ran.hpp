#pragma once

#include <cmath>

#include "doubles.hpp"

struct Ran {
  unsigned long long u, v, w;

  Ran(unsigned long long seed) { set_seed(seed); }

  unsigned long long int64() {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17;
    v ^= v << 31;
    v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);
    unsigned long long x = u ^ (u << 21);
    x ^= x >> 35;
    x ^= x << 4;
    return (x + v) ^ w;
  }

  void set_seed(unsigned long long seed) {
    v = 4101842887655102017LL;
    w = 1;
    u = seed ^ v;
    int64();
    v = u;
    int64();
    w = v;
    int64();
  }

  double U() { return 5.42101086242752217E-20 * (double)int64(); }

  double U(double min, double max) { return min + (max - min) * U(); }

  double N() {
    double uu, vv, xx, yy, qq;
    do {
      uu = U();
      vv = 1.7156 * (U() - 0.5);
      xx = uu - 0.449871;
      yy = std::abs(vv) + 0.386595;
      qq = xx * xx + yy * (0.19600 * yy - 0.25472 * xx);
    } while (qq > 0.27597 && (qq > 0.27846 || vv * vv > -4. * std::log(uu) * uu * uu));
    return vv / uu;
  }

  double N(double mean, double std) { return mean + std * N(); }

  double3 point_on_unit_sphere() {
#if 1
    double x1, x2, r2;
    do {
      x1 = U(-1.0, 1.0);
      x2 = U(-1.0, 1.0);
      r2 = x1 * x1 + x2 * x2;
    } while (r2 >= 1.0);
    double t = 2.0 * std::sqrt(1 - r2);
    return double3{x1 * t, x2 * t, 1.0 - 2.0 * r2};
#else
    double x = N();
    double y = N();
    double z = N();

    double r = std::sqrt(x * x + y * y + z * z);
    return double3{x / r, y / r, z / r};
#endif
  }

  void point_on_unit_circle(double &x, double &y) {
    double x1, x2, r2;
    do {
      x1 = U(-1.0, 1.0);
      x2 = U(-1.0, 1.0);
      r2 = x1 * x1 + x2 * x2;
    } while (r2 >= 1.0 || r2 == 0.0);
    x = (x1 * x1 - x2 * x2) / r2;
    y = 2 * x1 * x2 / r2;
  }

  double3 unit_perp_to(double3 a) {
    double3 n = to_unit(a);

    double3 u;
    if (std::abs(n.x) > std::abs(n.z)) {
      double t = 1.0 / std::sqrt(n.x * n.x + n.y * n.y);
      u = {-n.y * t, n.x * t, 0.0};
    } else {
      double t = 1.0 / std::sqrt(n.y * n.y + n.z * n.z);
      u = {0.0, -n.z * t, n.y * t};
    }

    double3 v = cross(n, u);
    double c, s;
    point_on_unit_circle(c, s);
    return c * u + s * v;
  }
};
