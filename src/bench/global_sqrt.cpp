#include "od_global.hpp"
// Experiment 4: merge orbit panels, retain all physical resonance boundaries.
int main(int argc, char **argv) {
  double scale = std::getenv("OD_SCALE") ? std::atof(std::getenv("OD_SCALE")) : 1;
  return od::run(
      "global_sqrt",
      [=](auto const &r, auto const &f, auto &c) {
        od::Settings s;
        s.cap *= scale;
        s.orbit_tol = 1e-10 * scale;
        s.tol = 1e-8 * scale;
        s.scan = 1;
        return od::global_integrate(r, f, s, c);
      },
      argc, argv);
}
