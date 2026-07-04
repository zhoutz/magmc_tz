#pragma once

#include <cstdio>
#include <utility>

struct UniformHunt {
  double x_min, x_max;
  int n;

  UniformHunt(double x_min, double x_max, int n) : x_min(x_min), x_max(x_max), n(n) {}

  auto operator()(double &x) {
    int i;
    double a;

    if (x <= x_min) {
      if (x < x_min)
        printf("UniformHunt::operator(): x = %g is out of bounds [%g, %g]: %d\n", x, x_min, x_max,
               n);
      i = 0;
      a = 0;
      x = x_min;
      return std::make_pair(i, a);
    }
    if (x >= x_max) {
      if (x > x_max)
        printf("UniformHunt::operator(): x = %g is out of bounds [%g, %g]: %d\n", x, x_min, x_max,
               n);
      i = n - 2;
      a = 1;
      x = x_max;
      return std::make_pair(i, a);
    }
    double t = (x - x_min) / (x_max - x_min) * (n - 1);
    i = int(t);
    a = t - i;
    return std::make_pair(i, a);
  }
};
