#include "od_common.hpp"
// ref.cpp algorithm extended to nonradial or arbitrary distributions:
// same DOPR5 orbit + dense output + small (0.01r) steps + raw spatial QAGS,
// same solve_quadratic opacity. Replace radial beta=0 event with local mu,
// include all supplied support edges, and continue after D=0 (possible re-entry).
int main(int argc, char **argv) {
  double scale = std::getenv("OD_REF_SCALE") ? std::atof(std::getenv("OD_REF_SCALE")) : 1;
  return od::run(
      "ref_general",
      [=](auto const &r, auto const &f, auto &c) {
        od::Settings s;
        s.cap = .01 * scale;
        s.tol = 1e-6 * scale;
        s.orbit_tol = 1e-10 * scale;
        s.scan = 1;
        s.reference = true;
        return od::integrate(r, f, s, c);
      },
      argc, argv);
}
