#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "hunt.hpp"

struct BField {
  double Bpole, R_NS;
  double Delta_phi, p, A, C;
  double mu_min, mu_max;
  int mu_num;
  std::vector<double> f, fp;

  BField(std::string fname, double Bpole_, double R_NS_) {
    Bpole = Bpole_;
    R_NS = R_NS_;
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

  auto getB(double r, double mu) {
    double mu_sign = (mu >= 0) ? 1 : -1;
    mu = std::abs(mu);

    UniformHunt mu_hunt(mu_min, mu_max, mu_num);
    auto [i, a] = mu_hunt(mu);

    double fval = (1 - a) * f[i] + a * f[i + 1];
    double fpval = (1 - a) * fp[i] + a * fp[i + 1];
    fpval *= mu_sign;
    double sth = std::sqrt(1 - mu * mu);

    double Br = -fpval;
    double Bth = p * fval / sth;
    double Bph = A * std::pow(fval, 1 / p) * Bth;

    double scale = 0.5 * Bpole * std::pow(R_NS / r, 2 + p);

    return std::make_tuple(scale * Br, scale * Bth, scale * Bph);
  }
};
