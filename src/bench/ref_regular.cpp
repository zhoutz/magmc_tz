#include "od_common.hpp"
// Fine-step reference after removing the same known square-root singularity.
// Compare against raw ref_general too; this is a convergence cross-check, not
// an independent proof of the shared opacity formula.
int main(int argc, char **argv) {
  double scale = std::getenv("OD_REF_SCALE") ? std::atof(std::getenv("OD_REF_SCALE")) : 1;
  return od::run(
      "ref_regular",
      [=](auto const &r, auto const &f, auto &c) {
        od::Settings s;
        s.cap = .01 * scale;
        s.tol = 1e-8 * scale;
        s.orbit_tol = 1e-11 * scale;
        s.scan = 1;
        s.reference = true;
        s.transformed = true;
        return od::integrate(r, f, s, c);
      },
      argc, argv);
}
