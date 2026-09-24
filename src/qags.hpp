#ifndef QAGS_HPP_INCLUDED
#define QAGS_HPP_INCLUDED

// C++23, standard-library-only, finite-interval adaptive quadrature.
// Algorithmic references (this is an independent implementation, not a GSL port):
//   QUADPACK / GSL QAGS: adaptive quadrature with extrapolation;
//   P. Wynn, Mathematical Tables and Other Aids to Computation 10 (1956), 91-96.
// GSL's finite-interval QAGS uses G10/K21; this implementation uses G7/K15.
//
// Usage:
//   auto r = qags::integrate([](double x) { return 1 / std::sqrt(x); },
//                            0.0, 1.0);
//   if (r.converged()) { /* r.value, r.absolute_error */ }
//
// Tolerance: absolute_error <= max(absolute_tolerance,
//                                relative_tolerance * abs(value)).
// Error estimates are heuristic, not rigorous bounds. Narrow unsampled features
// and divergent integrals cannot in general be detected by finite sampling.
// Endpoints are never evaluated. Split at known interior singularities before
// calling integrate (especially if a singularity coincides with a rule node).
// Infinite bounds are not accepted; apply a change of variables first.
// Callable exceptions and allocation failures propagate to the caller.
// All state is local to each call; no dependency on GSL or global workspace.

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace qags {

enum class status {
    success,
    invalid_input,
    max_subdivisions,
    roundoff_detected,
    interval_too_small,
    non_finite_integrand,
    non_finite_arithmetic,
    probable_divergence
};

[[nodiscard]] constexpr std::string_view status_message(status code) noexcept {
    switch (code) {
    case status::success: return "requested tolerance reached";
    case status::invalid_input: return "invalid bounds, tolerances, or interval limit";
    case status::max_subdivisions: return "maximum number of active subintervals reached";
    case status::roundoff_detected: return "roundoff prevents further improvement";
    case status::interval_too_small: return "subinterval has insufficient floating-point resolution";
    case status::non_finite_integrand: return "integrand returned NaN or infinity";
    case status::non_finite_arithmetic: return "quadrature arithmetic overflowed";
    case status::probable_divergence: return "extrapolation suggests divergence or very slow convergence";
    }
    return "unknown status";
}

template<std::floating_point Real = double>
struct options {
    Real absolute_tolerance = std::max(Real(1e-10L), Real(100) * std::numeric_limits<Real>::epsilon());
    Real relative_tolerance = std::max(Real(1e-8L), Real(100) * std::numeric_limits<Real>::epsilon());
    std::size_t max_intervals = 1000;
    bool use_extrapolation = true;
};

template<std::floating_point Real = double>
struct result {
    Real value = std::numeric_limits<Real>::quiet_NaN();
    Real absolute_error = std::numeric_limits<Real>::infinity();
    status code = status::invalid_input;
    std::size_t evaluations = 0;
    std::size_t intervals = 0;       // Number of accepted active leaves.
    std::size_t extrapolations = 0; // Number of epsilon-table evaluations.
    bool extrapolated = false;     // Whether the returned value uses epsilon.

    [[nodiscard]] constexpr bool converged() const noexcept {
        return code == status::success;
    }
    [[nodiscard]] constexpr explicit operator bool() const noexcept {
        return converged();
    }
};

namespace detail {

template<class Real>
inline constexpr Real machine_epsilon = std::numeric_limits<Real>::epsilon();

// Neumaier summation also supports replacing a leaf by adding its negative.
template<class Real>
struct compensated_sum {
    Real main = 0;
    Real correction = 0;
    void add(Real x) noexcept {
        const Real t = main + x;
        correction += std::abs(main) >= std::abs(x)
                    ? (main - t) + x : (x - t) + main;
        main = t;
    }
    [[nodiscard]] Real value() const noexcept { return main + correction; }
};

// Avoid intermediate overflow/underflow in x*y*z when the final product fits.
template<class Real>
[[nodiscard]] Real product(Real x, Real y, Real z) noexcept {
    if (x == 0 || y == 0 || z == 0) return Real(0);
    int ex = 0, ey = 0, ez = 0;
    const Real mx = std::frexp(x, &ex);
    const Real my = std::frexp(y, &ey);
    const Real mz = std::frexp(z, &ez);
    return std::scalbn((mx * my) * mz, ex + ey + ez);
}

template<class Real>
struct rule_result {
    Real value = 0, error = 0, absolute = 0, deviation = 0;
    status code = status::success;
};

template<class Real, class F>
[[nodiscard]] rule_result<Real> kronrod15(F& f, Real a, Real b,
                                         std::size_t& evaluations) {
    // Positive Kronrod nodes, followed by the origin. Odd indices belong to G7.
    static constexpr std::array<Real, 8> nodes{
        Real(0.991455371120812639206854697526329L),
        Real(0.949107912342758524526189684047851L),
        Real(0.864864423359769072789712788640926L),
        Real(0.741531185599394439863864773280788L),
        Real(0.586087235467691130294144838258730L),
        Real(0.405845151377397166906606412076961L),
        Real(0.207784955007898467600689403773245L), Real(0)};
    static constexpr std::array<Real, 8> kw{
        Real(0.022935322010529224963732008058970L),
        Real(0.063092092629978553290700663189204L),
        Real(0.104790010322250183839876322541518L),
        Real(0.140653259715525918745189590510238L),
        Real(0.169004726639267902826583426598550L),
        Real(0.190350578064785409913256402421014L),
        Real(0.204432940075298892414161999234649L),
        Real(0.209482141084727828012999174891714L)};
    static constexpr std::array<Real, 4> gw{
        Real(0.129484966168869693270611432679082L),
        Real(0.279705391489276667901467771423780L),
        Real(0.381830050505118944950369775488975L),
        Real(0.417959183673469387755102040816327L)};

    rule_result<Real> out;
    const Real center = std::midpoint(a, b);
    const Real width = b - a;
    const Real half = std::isfinite(width) ? width / 2 : b / 2 - a / 2;
    std::array<Real, 15> x{}, values{};
    x[14] = center;
    for (std::size_t j = 0; j < 7; ++j) {
        x[2 * j] = std::lerp(a, b, (Real(1) - nodes[j]) / 2);
        x[2 * j + 1] = std::lerp(a, b, (Real(1) + nodes[j]) / 2);
    }
    // Check all abscissae before calling f, so rounding never samples an endpoint.
    for (Real t : x) {
        if (!(a < t && t < b)) {
            out.code = status::interval_too_small;
            return out;
        }
    }
    Real scale = 0;
    for (std::size_t j = 0; j < x.size(); ++j) {
        ++evaluations;
        values[j] = static_cast<Real>(std::invoke(f, x[j]));
        if (!std::isfinite(values[j])) {
            out.code = status::non_finite_integrand;
            return out;
        }
        scale = std::max(scale, std::abs(values[j]));
    }
    if (scale == 0) return out;
    for (Real& y : values) y /= scale;
    Real k = kw[7] * values[14], g = gw[3] * values[14];
    Real absolute = kw[7] * std::abs(values[14]);
    for (std::size_t j = 0; j < 7; ++j) {
        const Real pair = values[2 * j] + values[2 * j + 1];
        k += kw[j] * pair;
        absolute += kw[j] * (std::abs(values[2 * j]) + std::abs(values[2 * j + 1]));
        if (j % 2 == 1) g += gw[j / 2] * pair;
    }
    const Real mean = k / 2;
    Real deviation = kw[7] * std::abs(values[14] - mean);
    for (std::size_t j = 0; j < 7; ++j)
        deviation += kw[j] * (std::abs(values[2 * j] - mean)
                              + std::abs(values[2 * j + 1] - mean));
    Real error = std::abs(k - g);
    // QUADPACK-style variation rescaling and a roundoff floor.
    if (deviation != 0 && error != 0)
        error = deviation * std::pow(std::min(Real(1), Real(200) * (error / deviation)), Real(1.5));
    error = std::max(error, Real(50) * machine_epsilon<Real> * absolute);
    out.value = product(k, half, scale);
    out.absolute = product(absolute, half, scale);
    out.deviation = product(deviation, half, scale);
    out.error = product(error, half, scale);
    if (!std::isfinite(out.value) || !std::isfinite(out.absolute)
        || !std::isfinite(out.deviation) || !std::isfinite(out.error))
        out.code = status::non_finite_arithmetic;
    return out;
}

// A bounded, explicitly recomputed Wynn epsilon table. The recurrence is
// e[-1,n] = 0, e[0,n] = s[n],
// e[k+1,n] = e[k-1,n+1] + 1 / (e[k,n+1] - e[k,n]).
// Even columns approximate the limit; odd columns are auxiliary quantities.
// Singular cells invalidate only their descendants, preserving lower orders.
// With a fixed capacity, O(capacity^2) work and O(capacity) space per update.
template<class Real>
class epsilon_table {
    static constexpr std::size_t capacity = 32;
    std::array<Real, capacity> sequence_{};
    std::array<Real, 3> history_{};
    std::size_t size_ = 0, history_size_ = 0;

public:
    struct estimate {
        Real value = 0;
        Real error = std::numeric_limits<Real>::infinity();
        bool available = false;
    };

    [[nodiscard]] estimate append(Real s) noexcept {
        if (size_ == capacity) {
            std::move(sequence_.begin() + 1, sequence_.end(), sequence_.begin());
            --size_;
        }
        sequence_[size_++] = s;
        estimate answer;
        if (size_ < 3) return answer;
        std::array<Real, capacity> previous{}, current = sequence_, next{};
        Real best_score = std::numeric_limits<Real>::infinity();
        for (std::size_t order = 1; order < size_; ++order) {
            const std::size_t count = size_ - order;
            for (std::size_t i = 0; i < count; ++i) {
                const Real delta = current[i + 1] - current[i];
                const Real magnitude = std::max(std::abs(current[i]), std::abs(current[i + 1]));
                next[i] = std::numeric_limits<Real>::quiet_NaN();
                if (std::isfinite(delta) && std::isfinite(previous[i + 1])
                    && std::abs(delta) > Real(4) * machine_epsilon<Real> * magnitude) {
                    const Real candidate = previous[i + 1] + Real(1) / delta;
                    if (std::isfinite(candidate)) next[i] = candidate;
                }
            }
            if (order % 2 == 0) {
                const std::size_t i = count - 1; // Newest entry in this column.
                if (std::isfinite(next[i])) {
                    const Real score = std::abs(next[i] - previous[i + 1])
                                     + std::abs(next[i] - previous[i]);
                    if (score < best_score) {
                        best_score = score;
                        answer.value = next[i];
                        answer.available = true;
                    }
                }
            }
            previous = current;
            current = next;
        }
        if (!answer.available) {
            history_size_ = 0;
            return answer;
        }
        // Do not claim convergence from a single accelerated estimate.
        if (history_size_ == history_.size()) {
            answer.error = 0;
            for (Real old : history_) answer.error += std::abs(answer.value - old);
            answer.error = std::max(answer.error,
                Real(50) * machine_epsilon<Real> * std::abs(answer.value));
        }
        history_[0] = history_[1];
        history_[1] = history_[2];
        history_[2] = answer.value;
        if (history_size_ < history_.size()) ++history_size_;
        return answer;
    }
};

template<class Real>
struct interval {
    Real a, b;
    rule_result<Real> rule;
    std::size_t depth;
    bool active = true;
};

template<class Real>
struct heap_entry {
    Real error;
    std::size_t id;
    [[nodiscard]] bool operator<(const heap_entry& other) const noexcept {
        return error < other.error || (error == other.error && id > other.id);
    }
};

} // namespace detail

// Real is deduced from the bounds (float/double/long double). The callable is
// invoked by reference, so mutable and move-only callables are supported.
template<class F, std::floating_point Real = double>
    requires std::invocable<F&, Real>
          && std::convertible_to<std::invoke_result_t<F&, Real>, Real>
[[nodiscard]] result<Real> integrate(F&& f, Real a, Real b, options<Real> opts = {}) {
    using detail::machine_epsilon;
    result<Real> out;
    if (!std::isfinite(a) || !std::isfinite(b)
        || !std::isfinite(opts.absolute_tolerance) || opts.absolute_tolerance < 0
        || !std::isfinite(opts.relative_tolerance) || opts.relative_tolerance < 0
        || (opts.absolute_tolerance == 0 && opts.relative_tolerance < Real(50) * machine_epsilon<Real>)
        || opts.max_intervals == 0
        || opts.max_intervals > std::numeric_limits<std::size_t>::max() / 2)
        return out;
    if (a == b) {
        out.value = out.absolute_error = 0;
        out.code = status::success;
        return out;
    }
    const Real sign = a < b ? Real(1) : Real(-1);
    if (a > b) std::swap(a, b);
    const auto tolerance = [&](Real value) {
        return std::max(opts.absolute_tolerance, opts.relative_tolerance * std::abs(value));
    };
    auto first = detail::kronrod15<Real>(f, a, b, out.evaluations);
    if (first.code != status::success) {
        out.code = first.code;
        return out;
    }

    std::vector<detail::interval<Real>> leaves;
    leaves.reserve(std::min(opts.max_intervals, std::size_t(1024)));
    leaves.push_back({a, b, first, 0});
    std::priority_queue<detail::heap_entry<Real>> largest, eligible;
    largest.push({first.error, 0});
    eligible.push({first.error, 0});
    std::vector<std::size_t> deferred;
    std::size_t frontier = 1, stagnant_splits = 0, rejected_extrapolations = 0;
    detail::compensated_sum<Real> area, errors, absolute, large_errors;
    area.add(first.value);
    errors.add(first.error);
    absolute.add(first.absolute);
    large_errors.add(first.error);
    detail::epsilon_table<Real> epsilon;
    (void)epsilon.append(first.value);
    std::size_t sequence_size = 1;
    Real extrapolated_value = 0;
    Real extrapolated_error = std::numeric_limits<Real>::infinity();
    out.intervals = 1;

    const auto clean = [&](auto& heap) {
        while (!heap.empty() && !leaves[heap.top().id].active) heap.pop();
    };
    const auto finish = [&](status code, bool allow_extrapolation) {
        out.value = sign * area.value();
        out.absolute_error = std::max(Real(0), errors.value());
        if (allow_extrapolation && (extrapolated_error < out.absolute_error
            || (code == status::success && extrapolated_error <= tolerance(extrapolated_value)))) {
            out.value = sign * extrapolated_value;
            out.absolute_error = extrapolated_error;
            out.extrapolated = true;
        }
        out.code = code;
        return out;
    };

    for (;;) {
        const Real direct = area.value();
        const Real direct_error = std::max(Real(0), errors.value());
        const Real target = tolerance(direct);
        if (direct_error <= target && (out.intervals > 1 || first.error != first.deviation || first.error == 0))
            return finish(status::success, false);
        if (out.intervals == 1 && direct_error <= Real(100) * machine_epsilon<Real> * first.absolute)
            return finish(status::roundoff_detected, false);

        clean(largest);
        std::size_t selected = largest.top().id;
        if (opts.use_extrapolation && leaves[selected].depth >= frontier) {
            clean(eligible);
            // Complete the coarser part of this refinement epoch before using
            // its total integral as a sequence element. Both choices use heaps
            // ordered by error; only the extrapolation phase restricts depth.
            const Real coarse_error = std::max(Real(0), large_errors.value());
            if (!eligible.empty() && coarse_error > target) {
                selected = eligible.top().id;
            } else {
                const auto accelerated = epsilon.append(direct);
                if (++sequence_size >= 3) ++out.extrapolations;
                if (accelerated.available && std::isfinite(accelerated.error)) {
                    const Real floor = Real(50) * machine_epsilon<Real> * std::max(Real(0), absolute.value());
                    const Real error = std::max(accelerated.error, floor) + coarse_error;
                    // Epsilon can assign a finite analytic continuation to a
                    // divergent sequence. Reject gross disagreement unless the
                    // integral is dominated by cancellation.
                    const bool cancellation = std::max(std::abs(direct), std::abs(accelerated.value))
                                            < Real(0.01) * first.absolute;
                    const Real ratio = direct == 0 ? Real(1) : accelerated.value / direct;
                    const bool plausible = cancellation || (ratio > Real(0.01) && ratio < Real(100)
                                                            && direct_error <= std::abs(direct));
                    if (plausible && error < extrapolated_error) {
                        extrapolated_value = accelerated.value;
                        extrapolated_error = error;
                    }
                    if (!plausible) ++rejected_extrapolations;
                    if (plausible && error <= tolerance(accelerated.value)) {
                        extrapolated_value = accelerated.value;
                        extrapolated_error = error;
                        return finish(status::success, true);
                    }
                }
                ++frontier;
                for (std::size_t id : deferred)
                    if (leaves[id].active) eligible.push({leaves[id].rule.error, id});
                deferred.clear();
                large_errors = errors;
                continue;
            }
        }
        if (out.intervals >= opts.max_intervals)
            return finish(rejected_extrapolations >= 3 ? status::probable_divergence : status::max_subdivisions, true);

        const auto parent = leaves[selected]; // Copy before vector reallocation.
        const Real middle = std::midpoint(parent.a, parent.b);
        if (!(parent.a < middle && middle < parent.b))
            return finish(status::interval_too_small, true);
        const auto left = detail::kronrod15<Real>(f, parent.a, middle, out.evaluations);
        if (left.code != status::success) return finish(left.code, false);
        const auto right = detail::kronrod15<Real>(f, middle, parent.b, out.evaluations);
        if (right.code != status::success) return finish(right.code, false);
        // Commit only after both children have valid estimates. On failure the
        // returned value still covers the complete original interval.
        auto next_area = area, next_errors = errors, next_absolute = absolute;
        next_area.add(-parent.rule.value);
        next_area.add(left.value);
        next_area.add(right.value);
        next_errors.add(-parent.rule.error);
        next_errors.add(left.error);
        next_errors.add(right.error);
        next_absolute.add(-parent.rule.absolute);
        next_absolute.add(left.absolute);
        next_absolute.add(right.absolute);
        if (!std::isfinite(next_area.value()) || !std::isfinite(next_errors.value())
            || !std::isfinite(next_absolute.value()))
            return finish(status::non_finite_arithmetic, false);

        const Real child_error = left.error + right.error;
        const Real change = std::abs((left.value - parent.rule.value) + right.value);
        if (left.error != left.deviation && right.error != right.deviation
            && child_error >= Real(0.99) * parent.rule.error
            && change <= Real(32) * machine_epsilon<Real> * (left.absolute + right.absolute))
            ++stagnant_splits;

        area = next_area;
        errors = next_errors;
        absolute = next_absolute;
        leaves[selected].active = false;
        ++out.intervals;
        const std::size_t child_depth = parent.depth + 1;
        if (opts.use_extrapolation && parent.depth < frontier)
            large_errors.add(-parent.rule.error);
        for (const auto& child : {detail::interval<Real>{parent.a, middle, left, child_depth},
                                  detail::interval<Real>{middle, parent.b, right, child_depth}}) {
            const auto id = leaves.size();
            leaves.push_back(child);
            largest.push({child.rule.error, id});
            if (opts.use_extrapolation) {
                if (child_depth < frontier) {
                    eligible.push({child.rule.error, id});
                    large_errors.add(child.rule.error);
                } else {
                    deferred.push_back(id);
                }
            }
        }
        if (stagnant_splits >= 12 && errors.value() > tolerance(area.value()))
            return finish(status::roundoff_detected, true);
    }
}

} // namespace qags

#endif // QAGS_HPP_INCLUDED
