#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "doubles.hpp"
#include "hunt.hpp"

struct BField {
  double B_pole, R_star;
  double Delta_phi, p, A, C;
  double mu_min, mu_max;
  int mu_num;
  std::vector<double> f, fp;

  BField(std::string fname, double B_pole, double R_star) : B_pole(B_pole), R_star(R_star) {
    std::ifstream fin(fname);
    if (!fin) {
      std::cout << "Error: cannot open file " << fname << std::endl;
      std::exit(1);
    }
    fin >> Delta_phi >> p >> A >> C >> mu_min >> mu_max >> mu_num;
    f.resize(mu_num);
    fp.resize(mu_num);
    for (int i = 0; i < mu_num; i++) {
      fin >> f[i] >> fp[i];
    }
  }

  double3 calc_B(double r, double mu) const {
    double mu_sign = (mu >= 0) ? 1 : -1;
    mu = std::abs(mu);

    UniformHunt mu_hunt(mu_min, mu_max, mu_num);
    auto [i, a] = mu_hunt(mu);

    double fval = (1 - a) * f[i] + a * f[i + 1];
    double fpval = (1 - a) * fp[i] + a * fp[i + 1];
    fpval *= mu_sign;
    double sth = std::sqrt(1 - mu * mu);

    double Br = -fpval;
    double Bth = (sth == 0) ? 0 : p * fval / sth;
    double Bph = A * std::pow(fval, 1 / p) * Bth;

    double scale = 0.5 * B_pole * std::pow(R_star / r, 2 + p);

    return double3(scale * Br, scale * Bth, scale * Bph);
  }
};
