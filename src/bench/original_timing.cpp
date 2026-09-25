// Time the original algorithms without process startup or table-loading cost.
// No algorithm, tolerance, stopping rule or output is changed.
#define main archived_main
#ifdef BENCH_REFERENCE
#include "ref.cpp"
#else
#include "base.cpp"
#endif
#undef main
#include <chrono>
int main() {
#ifdef BENCH_REFERENCE
  const char *name = "ref_timed";
#else
  const char *name = "base_timed";
#endif
  auto start = std::chrono::steady_clock::now();
  FILE *fp = std::fopen((std::string("output/") + name + ".txt").c_str(), "w");
  if (!fp) return 1;
  for (double b : {-.1, -.2, -.3, -.4, -.5, -.6, -.7, -.8, -.9})
    for (double m : {0., .1, .2, .3, .4, .5, .6, .7, .8, .9})
      for (double e : {.01, .1, 1., 10., 100.})
        for (auto p : {Polarization::E, Polarization::O})
          std::fprintf(fp, "%.2f %.2f %.2f %d %.16e\n", b, m, e, p == Polarization::E,
                       total_optical_depth(b, m, e, p));
  std::fclose(fp);
  std::printf("%s radial cases=900 seconds=%.9f\n", name,
              std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
}
