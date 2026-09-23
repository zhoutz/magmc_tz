#pragma once

// Independent radial-ray reference for resonant.typ.  It integrates over the
// particle distribution, not d tau / d l, and never calls the production
// resonance-root solver, optical-depth kernel, RK stepper, f(beta), or Bessel
// normalization.  The same tabulated B field is deliberately used: the
// benchmark measures transport error, not magnetic-table/model error.
#include "../ode/bfield.hpp"
#include "../ode/constants.hpp"
#include "../ode/distribution.hpp"
#include "../ode/photon.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace resonance_reference {

struct Result {
  double tau = 0;
  double error = 0; // Quadrature + tail + roundoff estimate; excludes field-table error.
  std::size_t evaluations = 0;
};

namespace detail {
using Real = long double;

// Independent Gauss--Kronrod (7,15) quadrature.  Neither rule samples the
// endpoints, including distribution-support boundaries.
template <class F> std::pair<Real, Real> kronrod(F &f, Real a, Real b, std::size_t &evaluations) {
  static constexpr Real x[] = {
      0.991455371120812639206854697526329L, 0.949107912342758524526189684047851L,
      0.864864423359769072789712788640926L, 0.741531185599394439863864773280788L,
      0.586087235467691130294144838258730L, 0.405845151377397166906606412076961L,
      0.207784955007898467600689403773245L, 0};
  static constexpr Real wk[] = {
      0.022935322010529224963732008058970L, 0.063092092629978553290700663189204L,
      0.104790010322250183839876322541518L, 0.140653259715525918745189590510238L,
      0.169004726639267902826583426598550L, 0.190350578064785409913256402421014L,
      0.204432940075298892414161999234649L, 0.209482141084727828012999174891714L};
  static constexpr Real wg[] = {
      0.129484966168869693270611432679082L, 0.279705391489276667901467771423780L,
      0.381830050505118944950369775488975L, 0.417959183673469387755102040816327L};
  const Real middle = (a + b) / 2, half = (b - a) / 2;
  const Real fm = f(middle);
  Real k = wk[7] * fm, g = wg[3] * fm;
  for (int i = 0; i < 7; ++i) {
    const Real sum = f(middle - half * x[i]) + f(middle + half * x[i]);
    k += wk[i] * sum;
    if (i % 2 == 1) g += wg[i / 2] * sum;
  }
  evaluations += 15;
  return {half * k, std::abs(half * (k - g))};
}

template <class F>
std::pair<Real, Real> adaptive(F &f, Real a, Real b, Real atol, int depth,
                               std::size_t &evaluations) {
  auto [value, error] = kronrod(f, a, b, evaluations);
  if (error <= atol + 64 * std::numeric_limits<Real>::epsilon() * std::abs(value))
    return {value, error};
  if (depth == 0) throw std::runtime_error("radial reference quadrature did not converge");
  const Real middle = (a + b) / 2;
  auto left = adaptive(f, a, middle, atol / 2, depth - 1, evaluations);
  auto right = adaptive(f, middle, b, atol / 2, depth - 1, evaluations);
  return {left.first + right.first, left.second + right.second};
}
} // namespace detail

// Propagation from photon.r to r_end at constant colatitude. alpha must be 0
// (outgoing) or pi (incoming), with endpoints consistent with its sign.
//
// For q=2+p, L=sqrt(1-rs/r), Q=q-rs/[2(r-rs)], changing variables in the
// delta-function optical-depth integral gives
//
// tau = (p+1) pi (Bphi/Btheta) integral d beta f(beta)/|mean beta|
//                  * |e'|^2 * (1-beta*mu) / (L[r(beta)] Q[r(beta)]).
//
// r(beta) solves x(r)=gamma*(1-beta*mu). Each supported beta contributes once,
// so this expression has no singularity at merging beta roots. Introduce
// u=sqrt(a*(gamma-1)) and use the exact Boltzmann mean to cancel K1e:
// (f/|mean beta|)*|d beta/du| = 2*(a+u^2)/sqrt(2*a+u^2)*exp(-u^2).
// In particular, arbitrarily narrow low-temperature layers remain resolved.
inline Result radial(const BField &field, const Boltzmann &distribution, double rs,
                     const Photon &photon, double r_end, double tolerance = 1e-10) {
  using detail::Real;
  if (!(tolerance > 0) || !std::isfinite(tolerance) || !(rs >= 0) || !(photon.r > rs) ||
      !(r_end > rs) || !(photon.omega_inf > 0))
    throw std::invalid_argument("invalid radial reference domain");
  if (std::abs(std::sin(photon.alpha)) > 1e-12)
    throw std::invalid_argument("radial reference requires radial photon");
  const int radial_sign = std::cos(photon.alpha) >= 0 ? 1 : -1;
  if ((r_end - photon.r) * radial_sign < 0)
    throw std::invalid_argument("radial endpoint opposes photon direction");
  if (photon.r == r_end) return {};

  const Real rlo = std::min(photon.r, r_end), rhi = std::max(photon.r, r_end);
  const Real q = 2 + Real(field.p), rs_ld = rs;
  // This holds everywhere outside the stellar surface in the project. It also
  // makes the inversion unique and exposes unsupported near-horizon inputs.
  if (!(q - rs_ld / (2 * (rlo - rs_ld)) > 0))
    throw std::invalid_argument("radial reference requires monotonic x(r)");
  const double3 rhat = std::cos(photon.psi) * photon.e1 + std::sin(photon.psi) * photon.e2;
  const double muz = std::clamp(rhat.z, -1.0, 1.0);
  const double3 B = field.calc_B(field.R_star, muz);
  const Real mu = radial_sign * Real(B.x / B.length());
  const Real prefactor =
      (1 + Real(field.p)) * std::acos(Real(-1)) * Real(field.Bphi_over_Btheta(muz));
  if (prefactor == 0) return {};
  const Real log_amplitude = std::log(Real(B_to_omega) * Real(B.length()) / Real(photon.omega_inf));
  const Real log_R = std::log(Real(field.R_star));
  const Real log_rlo = std::log(rlo), log_rhi = std::log(rhi);
  auto log_x = [&](Real logr) {
    const Real r = std::exp(logr);
    return log_amplitude - q * (logr - log_R) + std::log1p(-rs_ld / r) / 2;
  };
  const Real xmin = log_x(log_rhi), xmax = log_x(log_rlo);
  const Real b0 = distribution.b0;
  const Real s0 = std::sqrt((1 - b0) * (1 + b0));
  const Real a = s0 * (1 + s0) / (b0 * b0); // independently computed temperature
  const Real sign = b0 > 0 ? 1 : -1;
  auto beta = [&](Real u) { return sign * u * std::sqrt(2 * a + u * u) / (a + u * u); };
  auto log_g = [&](Real u) { return std::log1p(u * u / a) + std::log1p(-beta(u) * mu); };

  // exp(-64) suppresses the discarded thermal tail far below the requested
  // benchmark tolerances. The conservative absolute tail bound is returned
  // with the quadrature estimate below (for any positive temperature).
  constexpr Real umax = 8;
  std::vector<Real> turning{0, umax};
  if (sign * mu > 0 && std::abs(mu) < 1) {
    const Real smu = std::sqrt((1 - mu) * (1 + mu));
    const Real uturn = std::sqrt(a * mu * mu / (smu * (1 + smu)));
    if (uturn > 0 && uturn < umax) turning.push_back(uturn);
  }
  std::sort(turning.begin(), turning.end());
  std::vector<Real> breaks = turning;
  // Split at every beta for which the resonance lies at a propagation
  // endpoint. Without these splits even velocity quadrature can miss a thin
  // interval when the photon starts or stops inside the resonance layer.
  for (std::size_t i = 1; i < turning.size(); ++i) {
    for (const Real target : {xmin, xmax}) {
      Real lo = turning[i - 1], hi = turning[i];
      Real flo = log_g(lo) - target, fhi = log_g(hi) - target;
      if ((flo > 0) == (fhi > 0)) continue;
      for (int j = 0; j < 90; ++j) {
        const Real middle = (lo + hi) / 2;
        if (middle == lo || middle == hi) break;
        const Real fm = log_g(middle) - target;
        if ((flo > 0) == (fm > 0)) {
          lo = middle;
          flo = fm;
        } else
          hi = middle;
      }
      breaks.push_back((lo + hi) / 2);
    }
  }
  std::sort(breaks.begin(), breaks.end());
  breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
  auto integrand = [&](Real u) {
    const Real b = beta(u), target = log_g(u);
    Real lo = log_rlo, hi = log_rhi;
    for (int j = 0; j < 90; ++j) {
      const Real middle = (lo + hi) / 2;
      if (middle == lo || middle == hi) break;
      if (log_x(middle) > target)
        lo = middle;
      else
        hi = middle;
    }
    const Real r = std::exp((lo + hi) / 2), lapse = std::sqrt(1 - rs_ld / r);
    const Real Q = q - rs_ld / (2 * (r - rs_ld));
    const Real doppler = 1 - b * mu, mup = (mu - b) / doppler;
    const Real overlap = photon.pol == Polarization::E ? Real(0.5) : mup * mup / 2;
    const Real density = 2 * (a + u * u) / std::sqrt(2 * a + u * u) * std::exp(-u * u);
    return prefactor * density * overlap * doppler / (lapse * Q);
  };
  Result result;
  Real sum = 0, error = 0;
  for (std::size_t i = 1; i < breaks.size(); ++i) {
    const Real lo = breaks[i - 1], hi = breaks[i], gm = log_g((lo + hi) / 2);
    if (!(hi > lo) || gm < xmin || gm > xmax) continue;
    auto part = detail::adaptive(integrand, lo, hi, Real(tolerance) / breaks.size(), 24,
                                 result.evaluations);
    sum += part.first;
    error += part.second;
  }
  // density <= (sqrt(2a)+2u) exp(-u^2); overlap*doppler <= 1.
  const Real min_lapse = std::sqrt(1 - rs_ld / rlo);
  const Real min_Q = q - rs_ld / (2 * (rlo - rs_ld));
  const Real tail_bound = prefactor * (std::sqrt(2 * a) / (2 * umax) + 1) * std::exp(-umax * umax) /
                          (min_lapse * min_Q);
  result.tau = double(sum);
  const Real roundoff =
      (64 * std::numeric_limits<Real>::epsilon() + 2 * std::numeric_limits<double>::epsilon()) *
      std::abs(sum);
  result.error = double(error + tail_bound + roundoff);
  return result;
}

} // namespace resonance_reference
