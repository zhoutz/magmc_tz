// velocity_panels 的单文件阅读版。
// 建议阅读顺序：main -> total_optical_depth -> integrate -> ResonancePanels。
// 执行流程：三维轨道步进 -> 在步内定位共振边界 -> 分段积分 dτ/dl -> 累加 τ。
// 本文件不包含方法选择、其他积分方法或通用 benchmark 框架。
// 只复用原项目的基础组件；本方法的求积、边界搜索及端点处理都在此文件中。
//
// clang-format off
// 在仓库根目录编译、运行：
// g++-16 src/bench/only_velocity_panels.cpp -o build/only_velocity_panels -std=c++23 -O3 -I/opt/homebrew/include -L/opt/homebrew/lib -lgsl -lgslcblas
// build/only_velocity_panels
// python py/cmp_od.py only_velocity_panels
// 输出：output/only_velocity_panels.txt（与 base.cpp 相同的 900 例及列顺序）。
// clang-format on

#include "../bfield.hpp"
#include "../constants.hpp"
#include "../distribution.hpp"
#include "../dopr5.hpp"
#include "../photon.hpp"
#include "../qags.hpp"
#include "../roots.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <print>
#include <stdexcept>
#include <utility>
#include <vector>

namespace velocity_panels {
constexpr double R_star = 10;
constexpr double rs = 1.4 * schwarzschild_radius_of_sun_in_km;
constexpr double B_pole = 1e14;

// 1. 方法需要的电子分布接口。
// pdf 是归一化的数密度分布，mean 是带符号的平均速度。
// knots 列出支持边界、不连续点，以及足以分辨各窄峰的速度标记。
// 算法不假设 Boltzmann、单侧或单峰；有限采样不能发现未提供标记的任意窄峰。
struct Distribution {
  std::function<double(double)> pdf;
  double low, high, mean;
  std::vector<double> knots;

  double f(double beta) const { return beta > low && beta < high ? pdf(beta) : 0.; }
};

// 2. 轨道和局部散射物理量。
// State = {r, psi, alpha}；τ 不进入 DOPR5 的状态向量。

using State = std::array<double, 3>;

inline State geodesic(State const &y) {
  double r = y[0], alpha = y[2], L = std::sqrt(1 - rs / r);
  return {L * std::cos(alpha), std::sin(alpha) / r,
          -std::sin(alpha) / (r * L) * (1 - 1.5 * rs / r)};
}

// logx = ln(ωc/ω)，mu = k̂·B̂，prefactor 是 dτ/dl 的公共前因子。
struct ResonancePoint {
  double logx, mu, prefactor;
  double discriminant() const { return std::expm1(2 * logx) + mu * mu; }
  // F_beta=0 表示该位置与给定 beta 共振。这是定位窄层的事件函数。
  double surface(double beta) const {
    return logx - std::log1p(-beta * mu) + .5 * std::log1p(-beta * beta);
  }
};

// at() 计算当前点的磁场、角度与红移；opacity() 对两支共振速度求和。
struct Medium {
  BField const &field;
  Distribution const &distribution;
  Photon const &photon;

  ResonancePoint at(State const &y) const {
    auto [r, psi, alpha] = y;
    double3 rhat = std::cos(psi) * photon.e1 + std::sin(psi) * photon.e2;
    // 避免单位向量的舍入误差使磁极的 mu_z 越出磁场插值表。
    double muz = std::clamp(rhat.z, -1., 1.);
    auto B = field.calc_B(r, muz);
    double norm = B.length();
    auto b = B / norm;
    double rho = std::hypot(rhat.x, rhat.y);
    double3 theta = rho > 1e-14 ? double3{rhat.x * muz / rho, rhat.y * muz / rho, -rho}
                                : double3{std::copysign(1., muz), 0, 0};
    double3 phi = rho > 1e-14 ? double3{-rhat.y / rho, rhat.x / rho, 0} : double3{0, 1, 0};
    double mu = b.x * std::cos(alpha) +
                std::sin(alpha) * (b.y * dot(photon.n, phi) - b.z * dot(photon.n, theta));
    return {std::log(norm * B_to_omega / photon.omega_inf) + .5 * std::log1p(-rs / r),
            std::clamp(mu, -1., 1.),
            (field.p + 1) * pi * field.Bphi_over_Btheta(muz) / (r * std::abs(distribution.mean))};
  }

  std::array<double, 2> branches(ResonancePoint const &p, double D) const {
    if (!(D > 0)) return {NAN, NAN};
    double x = std::exp(p.logx), den = 1 + D, s = x * std::sqrt(D);
    // 先计算无相消的一根，再用韦达定理恢复另一根。
    double large = (p.mu + std::copysign(s, p.mu)) / den;
    // 两根都必须使用端点处理后的同一个 D，避免其中一根保留相消误差。
    double small = (p.mu * p.mu - D) / (den * large);
    if (!std::isfinite(small)) small = (p.mu - std::copysign(s, p.mu)) / den;
    return {std::min(large, small), std::max(large, small)};
  }

  // 用 |mu-beta| = sqrt(D)*(1-beta*mu)/x 避免合并点附近的相消。
  double opacity(ResonancePoint const &p, double D, double low = -1, double high = 1) const {
    if (!(D > 0) || p.prefactor == 0) return 0;
    double x = std::exp(p.logx), s = std::sqrt(D), sum = 0;
    // 共振条件下 mu'^2=D/x^2；对两个根分别取 PDF 并累加。
    double overlap = photon.pol == Polarization::E ? .5 : .5 * D / (x * x);
    for (double b : branches(p, D)) {
      if (b < low || b > high) continue;
      double f = distribution.f(b);
      if (f == 0) continue;
      sum += f * overlap * (1 - b * p.mu) * (1 - b * b) * x / s;
    }
    double value = p.prefactor * sum;
    if (!std::isfinite(value) || value < 0) throw std::runtime_error("Invalid resonant opacity");
    return value;
  }

  // 沿完整轨道求方向导数，包含磁场、光子方向和引力红移的变化。
  // 本方法只在 D=0 端点的 Taylor 展开中使用 rates，不用它限制轨道步长。
  std::array<double, 2> rates(State const &y) const {
    auto dy = geodesic(y);
    double h = 1e-5 * y[0];
    State a = y, b = y;
    for (int i = 0; i < 3; ++i) {
      a[i] -= h * dy[i];
      b[i] += h * dy[i];
    }
    auto pa = at(a), pb = at(b);
    return {(pb.logx - pa.logx) / (2 * h), (pb.mu - pa.mu) / (2 * h)};
  }
};

inline Photon make_photon(double r, double muz, double alpha, double azimuth, double energy,
                          Polarization pol) {
  double sth = std::sqrt(std::max(0., 1 - muz * muz));
  double3 er{sth, 0, muz}, etheta{muz, 0, -sth}, ephi{0, 1, 0};
  double3 tangent = std::cos(azimuth) * etheta + std::sin(azimuth) * ephi;
  return {cross(er, tangent), er, tangent, r, 0, alpha, energy, pol};
}

// 3. 每个已分割区间的数值求积：开节点 Gauss 8/16 点对。
// 节点不取区间端点，避免直接求值可积奇点；两阶差仅是经验误差估计。
// 注意：它不能自行发现未知窄峰，因此必须先做第 4 部分的事件分段。
namespace quadrature {

template <class F> double gauss(F const &f, double a, double b, bool high) {
  static constexpr double x8[] = {.18343464249564980494, .52553240991632898582,
                                  .79666647741362673959, .96028985649753623168};
  static constexpr double w8[] = {.36268378337836198297, .31370664587788728734,
                                  .22238103445337447054, .10122853629037625915};
  static constexpr double x16[] = {
      .095012509837637440185, .28160355077925891323, .45801677765722738634, .61787624440264374845,
      .75540440835500303390,  .86563120238783174388, .94457502307323257608, .98940093499164993260};
  static constexpr double w16[] = {.18945061045506849629,  .18260341504492358887,
                                   .16915651939500253819,  .14959598881657673208,
                                   .12462897125553387205,  .095158511682492784810,
                                   .062253523938647892863, .027152459411754094852};
  const double *x = high ? x16 : x8, *w = high ? w16 : w8;
  int n = high ? 8 : 4;
  double m = (a + b) / 2, h = (b - a) / 2, sum = 0;
  for (int i = 0; i < n; ++i) sum += w[i] * (f(m - h * x[i]) + f(m + h * x[i]));
  return h * sum;
}
template <class F>
double adaptive(F const &f, double a, double b, double atol, double rtol, int depth = 0) {
  if (a == b) return 0;
  double hi = gauss(f, a, b, true), lo = gauss(f, a, b, false);
  if (!std::isfinite(hi) || !std::isfinite(lo)) throw std::runtime_error("Nonfinite quadrature");
  if (std::abs(hi - lo) <= atol + rtol * std::abs(hi)) return hi;
  if (depth >= 40) throw std::runtime_error("Quadrature did not converge");
  double m = (a + b) / 2;
  return adaptive(f, a, m, atol / 2, rtol, depth + 1) +
         adaptive(f, m, b, atol / 2, rtol, depth + 1);
}
} // namespace quadrature

// 4. 一个轨道步内的事件搜索及分段积分。
// 路径参数 t∈[0,1]，实际路径长度 l=l_old+t*h；不以半径作为积分变量。
// singular=true 标记 D=0 的根，积分时需要端点正则化。
struct Boundary {
  double t;
  bool singular = false;
};

// A path segment is parameterized by t in [0,1], NOT by radius. Events
// remain valid through radial turning points and changing magnetic angles.
template <class Path> struct ResonancePanels {
  Medium const &medium;
  Path const &path;
  double length, low, high;
  std::vector<Boundary> boundaries;

  ResonancePoint point(double t) const { return medium.at(path(t)); }
  // beta=2 是内部标记（不可能是电子速度）：表示 D=0 的事件函数。
  // ln(x)-0.5*ln(1-mu^2) 与 D 的符号相同，避免大 x 时使用巨大判别式。
  double surface(double t, double beta) const {
    auto p = point(t);
    return beta == 2 ? p.logx - .5 * std::log(std::max(1e-30, 1 - p.mu * p.mu)) : p.surface(beta);
  }

  ResonancePanels(Medium const &m, Path const &p, double h, std::vector<double> const &velocities,
                  int probes, double lo = -1, double hi = 1)
      : medium(m), path(p), length(h), low(lo), high(hi) {
    // 4a. 先采样轨道步，再对每个速度标记寻找所有可见的交点。
    boundaries = {{0, false}, {1, false}};
    std::vector<ResonancePoint> samples;
    for (int i = 0; i <= probes; ++i) samples.push_back(point(double(i) / probes));
    std::vector<double> surfaces = velocities;
    surfaces.push_back(2); // D=0
    for (double beta : surfaces) {
      std::vector<std::pair<double, double>> nodes;
      for (int i = 0; i <= probes; ++i) {
        auto p = samples[i];
        nodes.push_back({double(i) / probes,
                         beta == 2 ? p.logx - .5 * std::log(std::max(1e-30, 1 - p.mu * p.mu))
                                   : p.surface(beta)});
      }
      // 端点同号不意味着内部没有窄层。若采样斜率变号，先用黄金分割
      // 搜索极值，再把极值作为新的括根节点。有限探测不是全局完备证明；
      // 换分布或轨道后仍需检查 probes 和步长的收敛性。
      // Endpoints can have equal signs around a thin resonance pocket.
      // Locate interior extrema before root bracketing, in addition to
      // sampling. Probe refinement is independently checked in validation.
      for (int i = 1; i < probes; ++i) {
        double s1 = nodes[i].second - nodes[i - 1].second,
               s2 = nodes[i + 1].second - nodes[i].second;
        if (s1 * s2 >= 0) continue;
        double a = nodes[i - 1].first, b = nodes[i + 1].first, sgn = s1 < 0 ? 1. : -1.;
        constexpr double g = .6180339887498948482;
        double c = b - g * (b - a), d = a + g * (b - a), fc = sgn * surface(c, beta),
               fd = sgn * surface(d, beta);
        for (int k = 0; k < 48; ++k) {
          if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - g * (b - a);
            fc = sgn * surface(c, beta);
          } else {
            a = c;
            c = d;
            fc = fd;
            d = a + g * (b - a);
            fd = sgn * surface(d, beta);
          }
        }
        double t = (a + b) / 2;
        nodes.push_back({t, surface(t, beta)});
      }
      std::sort(nodes.begin(), nodes.end());
      for (size_t i = 1; i < nodes.size(); ++i) {
        auto [a, fa] = nodes[i - 1];
        auto [b, fb] = nodes[i];
        if (a == b) continue;
        if (fa == 0) boundaries.push_back({a, beta == 2});
        if (fb == 0) boundaries.push_back({b, beta == 2});
        if (std::signbit(fa) == std::signbit(fb)) continue;
        double t = zriddr([&](double t) { return surface(t, beta); }, a, b, 3e-15);
        boundaries.push_back({t, beta == 2});
      }
    }
    // 4b. 把所有速度边界与 D=0 交点合并、排序、去重。
    std::sort(boundaries.begin(), boundaries.end(), [](auto a, auto b) { return a.t < b.t; });
    std::vector<Boundary> clean;
    for (auto b : boundaries) {
      if (!clean.empty() && b.t - clean.back().t < 8e-15)
        clean.back().singular |= b.singular;
      else
        clean.push_back(b);
    }
    boundaries = std::move(clean);
  }

  // 4c. 只在已经分好的区间 [left.t, right.t] 上积分 τ。
  // 找到边界后，才可用区间中点判断是否存在支持范围内的共振根。
  double integrate(Boundary left, Boundary right, double atol, double rtol) {
    double a = left.t, b = right.t, w = b - a;
    if (!(w > 0)) return 0;
    auto mid = point((a + b) / 2);
    double D = mid.discriminant();
    if (D <= 0) return 0; // all D=0 crossings have been split
    int contributing = 0;
    for (double beta : medium.branches(mid, D))
      contributing += beta > std::max(low, medium.distribution.low) &&
                      beta < std::min(high, medium.distribution.high);
    if (contributing == 0) return 0;

    // Local Taylor coefficients anchor D at a located simple caustic.
    // They avoid subtracting nearly equal O(1) numbers in tiny neighborhoods
    // after the square-root change of variable. No angular specialization.
    // D=0 端点附近直接计算 x^2+mu^2-1 会丢失有效位数。
    // 用沿完整轨道的 D 一、二阶导数，以已定位的根为零点作局部展开。
    // 一阶系数带 direction，二阶系数不带；再乘 h、h^2 转为 t 的导数。
    auto coefficients = [&](double t, double direction) {
      auto y = path(t), p = medium.at(y);
      auto rates = medium.rates(y);
      double first = 2 * std::exp(2 * p.logx) * rates[0] + 2 * p.mu * rates[1];
      double dl = 1e-4 * y[0];
      auto dy = geodesic(y);
      State yp = y, ym = y;
      for (int k = 0; k < 3; ++k) {
        yp[k] += dl * dy[k];
        ym[k] -= dl * dy[k];
      }
      auto pp = medium.at(yp), pm = medium.at(ym);
      auto rp = medium.rates(yp), rm = medium.rates(ym);
      double second = (2 * std::exp(2 * pp.logx) * rp[0] + 2 * pp.mu * rp[1] -
                       2 * std::exp(2 * pm.logx) * rm[0] - 2 * pm.mu * rm[1]) /
                      (2 * dl);
      return std::array<double, 2>{direction * first * length, second * length * length};
    };
    std::array<double, 2> lc{}, rc{};
    if (left.singular) lc = coefficients(a, 1);
    if (right.singular) rc = coefficients(b, -1);
    // t=a+(b-a)*sin^2(pi*z/2)，dt/dz=(b-a)*pi*sin*cos。
    // E 模的 D^(-1/2) 与 Jacobian 的零点相消；还必须乘实际步长 h。
    auto integrand = [&](double z) {
      // sin^2 maps both possible endpoint singularities at once.
      double sa = std::sin(.5 * pi * z), sb = std::cos(.5 * pi * z);
      double da = w * sa * sa, db = w * sb * sb;
      double t = z < .5 ? a + da : b - db;
      auto p = point(t);
      double D = p.discriminant();
      // A panel may be extremely narrow in D. Scaling the Taylor region by
      // panel width alone leaves most of such a panel dominated by roundoff.
      if (left.singular && da * length < 1e-5 * path(a)[0] && std::abs(lc[0] * da) < 1e-5)
        D = lc[0] * da + .5 * lc[1] * da * da;
      if (right.singular && db * length < 1e-5 * path(b)[0] && std::abs(rc[0] * db) < 1e-5)
        D = rc[0] * db + .5 * rc[1] * db * db;
      return medium.opacity(p, D, low, high) * length * w * pi * sa * sb;
    };
    return quadrature::adaptive(integrand, 0, 1, atol, rtol);
  }
};

// 5. 唯一的传输主循环（没有 Method 枚举，也没有其他方法的分支）。
struct Options {
  double tolerance = 1e-6;
  double resolution = .2;       // 几何步长上限 h/r，内部另有 0.4 的保护上限
  double initial_step = .01;    // km
  double escape_radius = 10000; // km
  int probes = 8;               // 每个轨道步的事件探测分段数
  double path_limit = INFINITY; // 可选：只求一段路径的累计深度
  double tau_target = INFINITY; // 可选：指定散射光学深度
};
struct Result {
  double tau = 0, length = 0;
  State state{};
  int outcome = 0; // 0=逃逸，1=撞星，2=达到路径长度上限，3=散射
};
Result integrate(BField const &field, Distribution const &distribution, Photon const &photon,
                 Options const &o = {}) {
  if (!(o.tolerance > 0 && o.resolution > 0 && o.initial_step > 0 && o.probes >= 4))
    throw std::invalid_argument("Invalid integration options");
  if (!(std::isfinite(distribution.mean) && distribution.mean != 0))
    throw std::invalid_argument("Current-normalized opacity requires a nonzero mean velocity");
  if (!(photon.r >= R_star && photon.alpha >= 0 && photon.alpha <= pi && photon.omega_inf > 0 &&
        o.path_limit > 0 && o.tau_target > 0 && o.escape_radius > photon.r))
    throw std::invalid_argument("Invalid photon or termination bounds");
  Result result;
  result.state = {photon.r, photon.psi, photon.alpha};
  if (photon.r == R_star && photon.alpha > pi / 2) {
    result.outcome = 1;
    return result;
  }
  Medium medium{field, distribution, photon};
  // 5a. 准备要追踪的速度边界。±1 不可代入 log(1-beta^2)。
  std::vector<double> velocities = distribution.knots;
  if (distribution.low > -1) velocities.push_back(distribution.low);
  if (distribution.high < 1) velocities.push_back(distribution.high);
  velocities.push_back(0);
  for (double v : velocities)
    if (!std::isfinite(v) || v < -1 || v > 1)
      throw std::invalid_argument("Velocity landmarks must lie in [-1,1]");
  std::erase_if(velocities, [](double v) { return std::abs(v) == 1; });
  std::sort(velocities.begin(), velocities.end());
  velocities.erase(std::unique(velocities.begin(), velocities.end()), velocities.end());
  // 5b. DOPR5 只积分三个轨道变量，不让平坦的 τ 决定轨道步长。
  auto rhs = [](double, State const &y, State &dy) { dy = geodesic(y); };
  double orbit_tol = std::min(1e-10, o.tolerance * .01);
  StepperDopr5<3, decltype(rhs)> s(rhs, orbit_tol, orbit_tol);
  s.add_event([&](double, State const &y) { return y[0] - o.escape_radius; });
  s.add_event([](double, State const &y) { return y[0] - R_star; });
  s.init(0, o.initial_step, result.state);
  long steps = 0;
  while (true) {
    if (++steps > 2000000) throw std::runtime_error("Orbit step budget exceeded");
    double r = s.y_old[0];
    double cap = o.resolution * r;
    cap = std::min({cap, .4 * r, o.path_limit - s.x_old});
    if (cap <= 0) {
      result.outcome = 2;
      break;
    }
    // 5c. 默认 h<=0.2r；防漏层主要依靠后面的步内搜索，而非此几何限步。
    s.do_step(cap);
    int event = s.detect_event();
    // detect_event 若截短了步长，会保留原步的 dense polynomial。
    // 此时不能再次 prepare_dense，否则插值的归一化步长会被改错。
    if (event < 0) s.prepare_dense();
    // 即使一个步的两端都在星外，也要检查中途的半径极小值是否已进入星体。
    if (s.y_old[2] > pi / 2 && s.y_new[2] <= pi / 2) {
      double turn = zriddr([&](double l) { return s.dense_out(l)[2] - pi / 2; }, s.x_old,
                           s.x_old + s.h_old, 1e-12);
      if (s.dense_out(turn)[0] <= R_star) {
        double hit =
            zriddr([&](double l) { return s.dense_out(l)[0] - R_star; }, s.x_old, turn, 1e-12);
        s.h_old = hit - s.x_old;
        s.y_new = s.dense_out(hit);
        event = 1;
      }
    }
    // 5d. 在整段轨道的 dense output 内查找速度事件，再逐区间累加 τ。
    double h = s.h_old;
    auto path = [&](double t) { return s.dense_out(s.x_old + t * h); };
    ResonancePanels panels(medium, path, h, velocities, o.probes, distribution.low,
                           distribution.high);
    for (size_t i = 1; i < panels.boundaries.size(); ++i) {
      auto left = panels.boundaries[i - 1], right = panels.boundaries[i];
      double atol = o.tolerance * .01 * std::min(1., h / r) / (panels.boundaries.size() - 1);
      auto integrate_panel = [&](Boundary end) {
        return panels.integrate(left, end, atol, o.tolerance);
      };
      double dtau = integrate_panel(right);
      // 可选：定位给定光学深度处的散射事件。默认目标为无穷，不走此分支。
      // 对区间内的累计积分再次求根，返回实际路径位置，而非整个步的终点。
      if (result.tau + dtau >= o.tau_target) {
        double needed = o.tau_target - result.tau;
        double t = zriddr(
            [&](double t) {
              return t == right.t ? dtau - needed : integrate_panel({t, false}) - needed;
            },
            left.t, right.t, 1e-11);
        result.tau = o.tau_target;
        result.length = s.x_old + t * h;
        result.state = path(t);
        result.outcome = 3;
        return result;
      }
      result.tau += dtau;
    }
    result.length = s.x_old + h;
    result.state = s.y_new;
    // 经过 D=0 不终止：非径向光线可能重新进入共振区。
    // 只在逃逸、撞星或到达指定路径长度时终止。
    if (event >= 0) {
      result.outcome = event;
      break;
    }
    if (result.length >= o.path_limit) {
      result.outcome = 2;
      break;
    }
    s.update_old();
  }
  return result;
}

// 6. 本例的输入准备：仅此函数了解 Boltzmann 分布。
// QAGS 在这里仅用于计算分位点；沿光线路径的 τ 积分使用上面的 Gauss 8/16。
// 分位点和尾部标记只决定分段位置，不裁剪 PDF 的尾部。
Distribution make_boltzmann_distribution(double b0) {
  Boltzmann fb(b0), positive(std::abs(b0));
  auto quantile = [&](double probability) {
    return zriddr(
        [&](double beta) {
          return qags([&](double v) { return positive.f(v); }, 0., beta, 1e-13, 1e-12) -
                 probability;
        },
        0., 1., 1e-14);
  };
  std::vector<double> knots = {0};
  double tail = quantile(.999);
  for (double probability : {.001, .01, .1, .25, .5, .75, .9, .99, .999})
    knots.push_back(std::copysign(quantile(probability), b0));
  for (double fraction : {.1, .25, .5, .75, .9, .99, .999})
    knots.push_back(std::copysign(tail + (1 - tail) * fraction, b0));
  return {[fb](double beta) { return fb.f(beta); }, fb.b_min, fb.b_max, fb.b_bar(), knots};
}

// 保留与 base.cpp/ref.cpp 类似的入口；photon 本身可以是非径向光子。
// 要换电子分布，只需构造 Distribution{pdf, low, high, mean, knots}。
double total_optical_depth(BField const &field, Distribution const &distribution,
                           Photon const &photon) {
  return integrate(field, distribution, photon).tau;
}
} // namespace velocity_panels

// 7. 与 base.cpp/ref.cpp 相同的径向测试主程序。
// make_photon 的 alpha、azimuth 可直接改成非径向初始条件；核心算法不变。
int main() {
  using namespace velocity_panels;
  try {
    BField field("table/bfield_t10.txt", B_pole, R_star);
    std::FILE *fp = std::fopen("output/only_velocity_panels.txt", "w");
    if (!fp) throw std::runtime_error("Cannot open output/only_velocity_panels.txt");
    try {
      for (double b0 : {-.1, -.2, -.3, -.4, -.5, -.6, -.7, -.8, -.9}) {
        // 每种分布只准备一次速度标记，不在每个轨道步重新计算分位点。
        auto distribution = make_boltzmann_distribution(b0);
        for (double muz : {0., .1, .2, .3, .4, .5, .6, .7, .8, .9})
          for (double omega_inf : {.01, .1, 1., 10., 100.})
            for (Polarization pol : {Polarization::E, Polarization::O}) {
              auto photon = make_photon(R_star, muz, 0., 0., omega_inf, pol);
              double tau = total_optical_depth(field, distribution, photon);
              std::println(fp, "{:.2f} {:.2f} {:.2f} {} {:.16e}", b0, muz, omega_inf,
                           pol == Polarization::E ? 1 : 0, tau);
            }
      }
    } catch (...) {
      std::fclose(fp);
      throw;
    }
    if (std::fclose(fp) != 0) throw std::runtime_error("Failed to write optical-depth output");
    return 0;
  } catch (std::exception const &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
