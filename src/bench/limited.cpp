#include "od_common.hpp"
// Experiment 1: minimal base.cpp repair, conservative spatial cap.
int main(int argc, char **argv) {
  return od::run(
      "limited",
      [](auto const &r, auto const &f, auto &c) { return od::limited(r, f, .01, 1e-8, c); }, argc,
      argv);
}
