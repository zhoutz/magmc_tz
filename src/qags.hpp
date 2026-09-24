#pragma once

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <memory>
#include <type_traits>

// Finite-interval integration. Returns the integral, or aborts on failure.
// f must not throw. No copy of f is made; capturing/mutable lambdas work.
// Each invocation owns its workspace, including nested/parallel invocations.
// This wrapper does not modify GSL's process-wide error handler.
template <class F>
  requires std::is_invocable_r_v<double, F &, double>
[[nodiscard]] double qags(F &&f, double a, double b, double epsabs = 1e-10, double epsrel = 1e-8,
                          std::size_t limit = 1000, double *abserr = nullptr) {
  if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(epsabs) || !std::isfinite(epsrel) ||
      epsabs < 0.0 || epsrel < 0.0 || limit == 0) {
    std::fputs("qags: invalid bounds, tolerances, or limit\n", stderr);
    std::abort();
  }

  using Workspace =
      std::unique_ptr<gsl_integration_workspace, decltype(&gsl_integration_workspace_free)>;
  Workspace workspace(gsl_integration_workspace_alloc(limit), &gsl_integration_workspace_free);
  if (!workspace) {
    std::fputs("qags: workspace allocation failed\n", stderr);
    std::abort();
  }

  // An object holding a reference also handles const callables and functions.
  auto callable = [&f](double x) noexcept -> double { return std::invoke(f, x); };

  gsl_function integrand{};
  integrand.function = +[](double x, void *params) noexcept -> double {
    return (*static_cast<decltype(callable) *>(params))(x);
  };
  integrand.params = std::addressof(callable);

  double result = 0.0;
  double error = 0.0;
  const int status = gsl_integration_qags(&integrand, a, b, epsabs, epsrel, limit, workspace.get(),
                                          &result, &error);

  // Also check the status when the application disables GSL's default abort.
  if (status != GSL_SUCCESS || !std::isfinite(result) || !std::isfinite(error)) {
    std::fprintf(stderr, "qags failed: %s (status=%d, result=%.17g, abserr=%.17g)\n",
                 status == GSL_SUCCESS ? "non-finite output" : gsl_strerror(status), status, result,
                 error);
    std::abort();
  }

  if (abserr) *abserr = error;
  return result;
}
