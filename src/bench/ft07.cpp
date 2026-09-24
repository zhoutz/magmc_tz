// FT07, sections 3.3--3.4, equations (31), (38), (39).
// See ft07_cap in bench_transport.hpp for the velocity-space step controller.
#include "bench_driver.hpp"
int main(int argc, char **argv) {
  return benchmark_main(argc, argv, "ft07", transport::Method::ft07);
}
