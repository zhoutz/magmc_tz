#include "od_common.hpp"
// Experiment 3: event splitting plus desingularized endpoint quadrature.
int main(int argc, char **argv) {
  return od::run(
      "event_sqrt",
      [](auto const &r, auto const &f, auto &c) {
        od::Settings s;
        s.transformed = true;
        return od::integrate(r, f, s, c);
      },
      argc, argv);
}
