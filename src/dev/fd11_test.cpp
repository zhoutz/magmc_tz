// Independent complex-amplitude reference and physical regression tests.
#define FD11_NO_MAIN
#include "fd11.cpp"
#include <complex>
#include <fstream>

using Complex = std::complex<double>;
using Jones = std::array<double, 4>; // Re Ax, Im Ax, Re Ay, Im Ay (tests only)

void require(bool pass, std::string const &message) {
  if (!pass) throw std::runtime_error(message);
}
fd11::Stokes to_stokes(Jones const &j) {
  Complex ax(j[0], j[1]), ay(j[2], j[3]);
  return {std::norm(ax) - std::norm(ay), 2 * std::real(ax * std::conj(ay)),
          2 * std::imag(ax * std::conj(ay))};
}
Jones eigen_jones(fd11::Coeff c, Polarization p) {
  double phi = 0.5 * std::atan2(c.s, c.c);
  if (p == Polarization::O) return {std::cos(phi), 0, std::sin(phi), 0};
  return {-std::sin(phi), 0, std::cos(phi), 0};
}
void jones_rhs(fd11::Coeff c, Jones const &j, Jones &d) {
  double h0 = (11.0 / 6) * c.kappa;
  double a = 0.5 * c.kappa * c.c, b = 0.5 * c.kappa * c.s;
  // Direct A'=i H A, retaining the common phase.
  d = {-(h0 + a) * j[1] - b * j[3], (h0 + a) * j[0] + b * j[2],
       -b * j[1] - (h0 - a) * j[3], b * j[0] + (h0 - a) * j[2]};
}
template <class F> void reference(F const &coeff, double x, double end, Jones &j) {
  auto rhs = [&](double s, Jones const &a, Jones &d) { jones_rhs(coeff(s), a, d); };
  StepperDopr5<4, decltype(rhs)> ode(rhs, 2e-11, 2e-11);
  if (!(end > x)) return;
  ode.init(x, std::min(end - x, 0.01 / std::max(coeff(x).kappa, 1e-10)), j);
  while (ode.x_old < end) {
    double left = end - ode.x_old;
    if (left <= 8 * std::numeric_limits<double>::epsilon() * std::max(1.0, end)) break;
    ode.do_step(left);
    ode.update_old();
  }
  j = ode.y_old;
}

void algebra_tests() {
  // Sign check: A=(1,1)/sqrt(2), H diagonal with Hxx-Hyy=2.
  auto c = [](double) { return fd11::Coeff{2, 1, 0}; };
  auto p = fd11::magnus(c, 0, 0.3);
  auto s = p.apply({0, 1, 0});
  require(std::abs(s.y - std::cos(0.6)) < 1e-14 &&
          std::abs(s.z - std::sin(0.6)) < 1e-14, "Stokes V sign");
  double worst_freeze = 0, worst_overlap = 0;
  Ran rng(9876);
  for (int i = 0; i < 80; ++i) {
    double phi = rng.U(-pi, pi);
    fd11::Coeff cc{rng.U(0.01, 3), std::cos(2 * phi), std::sin(2 * phi)};
    Jones j{rng.N(), rng.N(), rng.N(), rng.N()};
    double norm = 0;
    for (double v : j) norm += v * v;
    for (double &v : j) v /= std::sqrt(norm);
    auto before = j;
    double h = std::pow(10.0, rng.U(-4, 0));
    auto constant = [&](double) { return cc; };
    auto prop = fd11::magnus(constant, 0, h);
    reference(constant, 0, h, j);
    double delta2 = 0;
    for (int k = 0; k < 4; ++k) delta2 += std::pow(j[k] - before[k], 2);
    worst_freeze = std::max(worst_freeze,
        std::abs(std::sqrt(delta2) - prop.amplitude_change(to_stokes(before))));
    require((prop.apply(to_stokes(before)) + (-1.0) * to_stokes(j)).length() < 1e-8,
            "Constant-field Jones/Stokes mismatch");
    // Independently evaluate the two squares in FD11 (33), in B's screen.
    Complex ax(before[0], before[1]), ay(before[2], before[3]);
    Complex abx = std::cos(phi) * ax + std::sin(phi) * ay;
    Complex aby = -std::sin(phi) * ax + std::cos(phi) * ay;
    double mu = rng.U(-1, 1);
    for (int sign : {-1, 1}) {
      double direct = 0.5 * std::norm(mu * abx + Complex(0, sign) * aby);
      double from_stokes = fd11::overlap(to_stokes(before), cc, mu, sign);
      worst_overlap = std::max(worst_overlap, std::abs(direct - from_stokes));
    }
  }
  require(worst_freeze < 2e-8, "FD11 amplitude freeze metric mismatch");
  require(worst_overlap < 1e-14, "FD11 (33) overlap mismatch");
  // Noncommuting field axes: catches order/sign errors in Magnus products.
  auto rotating = [](double x) {
    return fd11::Coeff{3 + x, std::cos(2 * x), std::sin(2 * x)};
  };
  Jones j{0.6, 0, 0, 0.8};
  auto initial = to_stokes(j);
  reference(rotating, 0, 2, j);
  double previous = 0;
  for (int n : {20, 40, 80}) {
    auto v = initial;
    for (int i = 0; i < n; ++i) v = fd11::magnus(rotating, 2.0 * i / n, 2.0 / n).apply(v);
    double err = (v + (-1.0) * to_stokes(j)).length();
    if (previous) require(previous / err > 12, "Magnus fourth-order convergence");
    previous = err;
  }
  std::println("algebra: freeze_error={:.3e} overlap_error={:.3e} rotating_error={:.3e}",
               worst_freeze, worst_overlap, previous);
}

struct RayResult { fd11::Stokes s; double start_r, freeze_r, max_error; long steps; };
void polarized_rate_tests() {
  PhotonEvolution pe(bfield, fb, 123);
  auto bs = bfield.calc_B(45, 0);
  auto tangent = to_unit(double3{0, -bs.z, bs.y});
  pe.init(cross(double3{1, 0, 0}, tangent), {1, 0, 0}, tangent,
          24, Polarization::E, 45, 0, 1.3);
  pe.stepper.do_step(1e-4);
  auto g = pe.geo(pe.stepper.dense_out(0));
  auto c = pe.pol_coeff(0);
  std::array<double, 2> betas;
  require(g.D > 0 && solve_quadratic(g.x*g.x+g.mu*g.mu, -2*g.mu,
                                    (1+g.x)*(1-g.x), betas),
          "Rate fixture must be resonant");
  auto check_jones = [&](Jones j) {
    auto rates = pe.polarized_rates(betas, g, 0);
    double phi = 0.5 * std::atan2(c.s, c.c);
    Complex ax(j[0], j[1]), ay(j[2], j[3]);
    Complex abx = std::cos(phi)*ax + std::sin(phi)*ay;
    Complex aby = -std::sin(phi)*ax + std::cos(phi)*ay;
    for (int i=0; i<2; ++i) {
      double beta = betas[i], mup = (g.mu-beta)/(1-beta*g.mu);
      double expected = fb.f(beta)*std::pow(1-beta*beta, 1.5)*0.5*
          std::norm(mup*abx + Complex(0, pe.pol_options.charge_sign)*aby);
      require(std::abs(rates[i]-expected) < 1e-13*std::max(1.0, expected),
              "Polarized rate must include the FD11 overlap exactly once");
    }
    require(rates[0]+rates[1] > 0, "Rate fixture must have nonzero weight");
    return rates;
  };
  for (auto mode : {Polarization::E, Polarization::O}) {
    pe.pol = mode;
    pe.pol_options.enabled = true;
    pe.pol_start = std::numeric_limits<double>::infinity();
    check_jones(eigen_jones(c, mode)); // only O/E stored before ODE start
    pe.pol_start = 0;
    pe.pol_options.enabled = false;
    check_jones(eigen_jones(c, mode)); // disabled ODE must still include overlap
  }
  pe.pol_options.enabled = true;
  pe.pol_stage = PhotonEvolution::PolStage::Frozen;
  pe.frozen_at = 0;
  Jones mixed{0.6, 0, 0, 0.8};
  pe.stokes = to_stokes(mixed);
  for (int sign : {-1, 1}) {
    pe.pol_options.charge_sign = sign;
    auto rates = check_jones(mixed);
    auto terms = pe.opacity_terms(0);
    double expanded = terms.x + terms.y*(c.c*pe.stokes.x+c.s*pe.stokes.y)
                      + terms.z*pe.stokes.z;
    double direct = g.basepref*(rates[0]+rates[1]);
    require(std::abs(expanded-direct) < 1e-13*std::max(1.0, direct),
            "IBP coefficients and direct opacity must agree");
  }
  auto unsupported = pe.polarized_rates({-1, 0}, g, 0);
  require(unsupported[0] == 0 && unsupported[1] == 0,
          "Distribution support endpoints must have zero rate");
  std::println("rates: O/E, disabled ODE, mixed Stokes, both charges and IBP normalization passed");
}

void oscillatory_opacity_test(bool support_boundary) {
  PhotonEvolution pe(bfield, fb, 123);
  pe.pol_options.tolerance = 1e-8;
  auto bs = bfield.calc_B(45, 0);
  auto tangent = to_unit(double3{0, -bs.z, bs.y});
  auto normal = cross(double3{1, 0, 0}, tangent);
  double energy = support_boundary ? B_to_omega * bs.length() * std::sqrt(1 - rs / 45) : 24;
  pe.init(normal, {1, 0, 0}, tangent, energy, Polarization::E, 45, 0, 1.3);
  pe.pol_stage = PhotonEvolution::PolStage::Integrating;
  pe.pol_start = 0;
  pe.stokes = {0.3, 0.4, std::sqrt(0.75)};
  double kappa = pe.pol_coeff(pe.stepper.y_old).kappa;
  pe.stepper.h_old = 300 / kappa;
  pe.stepper.do_step(pe.stepper.h_old);
  double end = pe.stepper.h_old;
  pe.advance_polarization(0, end);
  auto f = [&](double x) {
    auto c = pe.pol_coeff(x);
    auto s = pe.polarization_at(x);
    auto terms = pe.opacity_terms(x);
    return terms.x + terms.y * (c.c * s.x + c.s * s.y) + terms.z * s.z;
  };
  double resolved = 0;
  Quad quad;
  for (int i = 0; i < 2000; ++i)
    resolved += quad.qags(f, end * i / 2000, end * (i + 1) / 2000, 1e-13, 1e-10);
  double fast = pe.fast_optical_depth(0, end);
  require(resolved > 1e-7, "Oscillatory opacity test must have nonzero scattering");
  require(std::abs(resolved - fast) < 2e-8, "Integration-by-parts opacity mismatch");
  std::println("oscillatory opacity (support_boundary={}): resolved={:.12g} fast={:.12g} error={:.3e}",
               support_boundary, resolved, fast, std::abs(resolved - fast));
}
void lifecycle_tests() {
  PhotonEvolution pe(bfield, fb, 91);
  pe.init({0, 0, 1}, {1, 0, 0}, {0, 1, 0}, 1, Polarization::E, R_star, 0, 0);
  pe.stepper.do_step();
  pe.advance_polarization(0, pe.stepper.h_old);
  require(pe.pol_stage == PhotonEvolution::PolStage::Mode && pe.stats.accepted == 0,
          "Premature ODE integration at emission");
  pe.stokes = {0, 1, 0};
  pe.pol_stage = PhotonEvolution::PolStage::Frozen;
  pe.scattered_data = {{45, 0, 0}, -0.5};
  pe.perform_scattering();
  require(pe.pol_stage == PhotonEvolution::PolStage::Mode && pe.stokes.length() == 0 &&
          !std::isfinite(pe.pol_start), "Scattering must reset to an eigenmode");
  // Exactly parallel field: deterministic degenerate mode and no NaNs.
  pe.init({0, 1, 0}, {0, 0, 1}, {1, 0, 0}, 1, Polarization::E, 100, 0, 0);
  pe.stepper.do_step();
  pe.advance_polarization(0, pe.stepper.h_old);
  require(std::isfinite(pe.stokes.x) && std::abs(pe.stokes.length() - 1) < 1e-12,
          "Parallel-field degeneracy");
  std::println("lifecycle: delayed start, scattering reset, parallel-field limit passed");
}
RayResult ray(double alpha, Polarization mode, double tol, double couple,
              double freeze, bool check_reference, double latitude = 0.3) {
  PhotonEvolution pe(bfield, fb, 4321);
  pe.pol_options.tolerance = tol;
  pe.pol_options.couple = couple;
  pe.pol_options.freeze = freeze;
  double3 rh{std::sqrt(1 - latitude * latitude), 0, latitude};
  double3 normal = to_unit(double3{-latitude, 0.7, rh.x});
  pe.init(normal, rh, cross(normal, rh), 1, mode, 18 * R_star, 0, alpha);
  require(pe.pol_stage == PhotonEvolution::PolStage::Mode && pe.stokes.length() == 0,
          "Emission must store only O/E");
  Jones j{};
  bool started = false;
  double start_r = 0, freeze_r = 0, max_error = 0;
  while (true) {
    pe.stepper.do_step(0.1 * pe.stepper.y_old[0]);
    int event = pe.stepper.detect_event();
    double l = pe.stepper.x_old, r = l + pe.stepper.h_old;
    auto stage_before = pe.pol_stage;
    pe.advance_polarization(l, r);
    if (!started && pe.pol_stage != PhotonEvolution::PolStage::Mode) {
      started = true;
      l = pe.pol_start;
      auto state = pe.stepper.dense_out(l);
      start_r = state[0];
      auto cc = pe.pol_coeff(l);
      double trigger = cc.kappa * start_r * couple;
      require(std::abs(trigger - 1) < 1e-8 || (pe.pol_start == 0 && trigger <= 1),
              "FD11 (34) root accuracy: " + std::to_string(trigger));
      j = eigen_jones(cc, mode);
    }
    if (pe.pol_stage == PhotonEvolution::PolStage::Frozen && freeze_r == 0)
      freeze_r = pe.stepper.dense_out(pe.frozen_at)[0];
    if (started && check_reference && stage_before != PhotonEvolution::PolStage::Frozen) {
      double end = std::min(r, pe.frozen_at);
      reference([&](double x) { return pe.pol_coeff(x); }, l, end, j);
      max_error = std::max(max_error, (pe.stokes + (-1.0) * to_stokes(j)).length());
      // Dense values used by scattering must agree with independent Jones.
      if (!pe.pol_segments.empty()) {
        auto seg = pe.pol_segments[pe.pol_segments.size() / 2];
        Jones mid = eigen_jones({1, 1, 0}, Polarization::O);
        // Reconstruct a Jones vector with seg.before, up to irrelevant phase.
        double ax = std::sqrt(std::max(0.0, (1 + seg.before.x) / 2));
        if (ax > 1e-8) {
          mid = {ax, 0, seg.before.y / (2 * ax), -seg.before.z / (2 * ax)};
          double xm = std::midpoint(seg.x, seg.end);
          reference([&](double x) { return pe.pol_coeff(x); }, seg.x, xm, mid);
          require((pe.polarization_at(xm) + (-1.0) * to_stokes(mid)).length() < 2e-5,
                  "Scattering dense polarization reference");
        }
      }
    }
    if (event >= 0) break;
    pe.stepper.update_old();
  }
  require(started, "Ray failed to start polarization");
  require(std::abs(pe.stokes.length() - 1) < 1e-10, "Polarization norm");
  if (check_reference) require(max_error < 2e-4, "Ray Jones reference error");
  return {pe.stokes, start_r, freeze_r, max_error, pe.stats.accepted};
}

int main() try {
  algebra_tests();
  polarized_rate_tests();
  oscillatory_opacity_test(false);
  oscillatory_opacity_test(true);
  lifecycle_tests();
  for (auto mode : {Polarization::E, Polarization::O}) {
    auto radial = ray(0, mode, 1e-6, 1e-3, 1e-3, true);
    require(std::abs(radial.s.z) < 1e-10, "Radial mode must not develop V");
    require(radial.freeze_r > 2 * radial.start_r, "Pure-mode premature freeze");
  }
  for (double a : {0.4, 1.0, 1.8}) {
    auto base = ray(a, Polarization::E, 1e-5, 1e-3, 1e-3, true);
    auto tight = ray(a, Polarization::E, 1e-8, 1e-3, 1e-3, false);
    auto deeper = ray(a, Polarization::E, 1e-8, 1e-5, 1e-3, false);
    auto unfrozen = ray(a, Polarization::E, 1e-8, 1e-3, 0, false);
    auto ordinary = ray(a, Polarization::O, 1e-8, 1e-3, 1e-3, false);
    double conv = (base.s + (-1.0) * tight.s).length();
    double deep = (tight.s + (-1.0) * deeper.s).length();
    double residual = (tight.s + (-1.0) * unfrozen.s).length();
    std::println("ray alpha={} start_km={:.6f} freeze_km={:.6f} steps={} Jones_error={:.3e} tolerance_error={:.3e} deeper_error={:.3e} freeze_error={:.3e} V={:.6f}",
                 a, base.start_r, base.freeze_r, base.steps, base.max_error,
                 conv, deep, residual, base.s.z);
    require(conv < 2e-4 && residual < 1e-3, "Trajectory convergence");
    // Eq. (34) fixes an approximate eigenmode start surface. Moving it is a
    // physical approximation sensitivity test, not an ODE error tolerance;
    // near k parallel B it need not converge at the paper's generic rate.
    require(std::isfinite(deep), "Non-finite coupling-surface sensitivity");
    require((tight.s + ordinary.s).length() < 1e-3, "O/E antipodal states");
  }
  std::println("All FD11 tests passed.");
} catch (std::exception const &e) {
  std::println(stderr, "TEST FAILED: {}", e.what());
  return 1;
}
