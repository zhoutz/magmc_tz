#include "od_common.hpp"
// Experiment 2: orbit and opacity separated; explicitly split resonance events.
int main(int argc, char **argv) {
  return od::run(
      "event_qags",
      [](auto const &r, auto const &f, auto &c) { return od::integrate(r, f, {}, c); }, argc, argv);
}
