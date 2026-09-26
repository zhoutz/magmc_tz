#pragma once

#include "../photon.hpp"
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <gsl/gsl_errno.h>

// Shared input/output only; each solver retains its historical numerical policy.
namespace optical_bench {
struct Input {
  int id = 0;
  double b0, muz, alpha, az, energy, radius;
  int pol;
  Polarization polarization() const { return pol ? Polarization::E : Polarization::O; }
};
struct Counts {
  size_t geometry = 0, density = 0, steps = 0, quadrature = 0, warnings = 0;
  double preparation_seconds = 0;
};
using Clock = std::chrono::steady_clock;

template <class Solve, class Prepare>
int run(int argc, char **argv, std::string const &name, Solve solve, Prepare prepare) {
  try {
    gsl_set_error_handler_off();
    bool const radial = argc == 1;
    if (!radial && (argc < 3 || argc > 4))
      throw std::invalid_argument("Usage: executable [input.dat output.dat [default|matched|cached|audit]]");
    std::string const mode = argc == 4 ? argv[3] : "default";
    std::vector<Input> inputs;
    if (radial) {
      for (int k = 1; k <= 9; ++k)
        for (int m = 0; m < 10; ++m)
          for (double e : {.01, .1, 1., 10., 100.})
            for (int p : {1, 0})
              inputs.push_back({int(inputs.size()), -k / 10., m / 10., 0, 0, e, 10, p});
      std::filesystem::create_directories("output");
    } else {
      std::ifstream in(argv[1]);
      if (!in) throw std::runtime_error("Cannot open benchmark input");
      Input c;
      while (in >> c.id >> c.b0 >> c.muz >> c.alpha >> c.az >> c.energy >> c.radius >> c.pol)
        inputs.push_back(c);
      if (!in.eof() || inputs.empty()) throw std::runtime_error("Malformed or empty benchmark input");
    }
    std::ofstream out(radial ? "output/" + name + ".txt" : argv[2]);
    if (!out) throw std::runtime_error("Cannot open benchmark output");
    out << std::setprecision(17);
    auto setup_start = Clock::now();
    prepare(inputs, mode);
    std::cerr << "preparation_seconds="
              << std::chrono::duration<double>(Clock::now() - setup_start).count() << '\n';
    int failures = 0;
    for (auto const &c : inputs) {
      Counts counts;
      double tau = std::numeric_limits<double>::quiet_NaN();
      std::string error;
      auto start = Clock::now();
      try {
        tau = solve(c, mode, counts);
        if (!std::isfinite(tau) || tau < 0) throw std::runtime_error("Invalid optical depth");
      } catch (std::exception const &e) {
        tau = std::numeric_limits<double>::quiet_NaN();
        error = e.what();
      }
      double seconds = std::chrono::duration<double>(Clock::now() - start).count();
      if (!error.empty()) {
        ++failures;
        std::cerr << "case " << c.id << ": " << error << '\n';
      }
      if (radial)
        out << c.b0 << ' ' << c.muz << ' ' << c.energy << ' ' << c.pol << ' ' << tau << '\n';
      else
        out << c.id << ' ' << tau << ' ' << seconds << ' ' << !error.empty() << ' '
            << counts.geometry << ' ' << counts.density << ' ' << counts.steps << ' '
            << counts.quadrature << ' ' << counts.warnings << ' ' << counts.preparation_seconds << '\n';
    }
    // Batch mode records failures per row, allowing the whole suite to finish.
    return radial && failures ? 2 : 0;
  } catch (std::exception const &e) {
    std::cerr << name << ": " << e.what() << '\n';
    return 1;
  }
}
template <class Solve>
int run(int argc, char **argv, std::string const &name, Solve solve) {
  return run(argc, argv, name, solve, [](auto const &, auto const &) {});
}
} // namespace optical_bench
