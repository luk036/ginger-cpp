/**
 * @file execution_policy.hpp
 * @brief Internal execution policies for parallelizable root-finding solvers
 *
 * Strategy + Template Method decomposition of the Bairstow / Aberth solver
 * families:
 *
 *   - Each algorithm (even/autocorr Bairstow, plain/autocorr Aberth) supplies a
 *     per-root "step" functor.  A step performs ONE Newton correction for root
 *     `idx` against the current values of the other roots and reports the
 *     residual tolerance of that root.
 *
 *   - Each execution mode (sequential / Jacobi-MT / atomic-decoupled) supplies
 *     a policy whose `run()` owns the iteration loop, chunk scheduling,
 *     synchronization and result aggregation.  The policy passes the step the
 *     accessors `get(idx)`/`set(value)` bound to its own data layout (plain
 *     vector, frozen snapshot, or `std::atomic` buffer) together with a
 *     neighbor index iterable (index order, or round-robin for the atomic
 *     variant).
 *
 * This removes the per-mode copy-paste of the solver loops while keeping the
 * public API (pbairstow_even_st/_mt/_atomic, aberth(_mt/_atomic), ...)
 * unchanged.
 */

#pragma once

#include <algorithm>               // for max, min
#include <atomic>                  // for atomic
#include <cmath>                   // for abs
#include <complex>                 // for complex
#include <cstddef>                 // for size_t
#include <future>                  // for future
#include <ginger/aberth.hpp>       // for horner_eval_c
#include <ginger/config.hpp>       // for Options
#include <ginger/robin.hpp>        // for fun::Robin
#include <ginger/rootfinding.hpp>  // for Vec2, delta_scalar, horner, suppress_old
#include <ginger/thread_pool.hpp>  // for ginger::get_thread_pool
#include <utility>                 // for pair
#include <vector>                  // for vector

namespace ginger::detail {

    /// @brief Iterable over [0, n) excluding `skip`, in ascending index order.
    class index_range_excluding {
      public:
        index_range_excluding(std::size_t n, std::size_t skip) : n_(n), skip_(skip) {}

        class iterator {
          public:
            iterator(std::size_t i, std::size_t skip) : i_(i), skip_(skip) {}
            auto operator*() const -> std::size_t { return i_; }
            auto operator!=(const iterator& other) const -> bool { return i_ != other.i_; }
            auto operator++() -> iterator& {
                ++i_;
                if (i_ == skip_) {
                    ++i_;
                }
                return *this;
            }

          private:
            std::size_t i_;
            std::size_t skip_;
        };

        auto begin() const -> iterator { return iterator{skip_ == 0 ? std::size_t{1} : 0, skip_}; }
        auto end() const -> iterator { return iterator{n_, skip_}; }

      private:
        std::size_t n_;
        std::size_t skip_;
    };

    /// @brief Derivative coefficients of a polynomial (Horner-usable for P'(z)).
    /// @param[in] coeffs Polynomial coefficients (highest degree first)
    /// @return vector of degree coefficients: c'[i] = (degree - i) * coeffs[i]
    inline auto derivative_coeffs(const std::vector<double>& coeffs) -> std::vector<double> {
        const auto degree = coeffs.size() - 1;
        auto coeffs1 = std::vector<double>(degree);
        for (auto idx = std::size_t{0}; idx < degree; ++idx) {
            coeffs1[idx] = static_cast<double>(degree - idx) * coeffs[idx];
        }
        return coeffs1;
    }

    /// @brief Roots of a quadratic factor x^2 - r x - q (real or complex pair).
    /// @param[in] vr Quadratic factor (r, q)
    /// @return The two roots
    inline auto roots_from_quadratic(const Vec2& vr)
        -> std::pair<std::complex<double>, std::complex<double>> {
        const auto r = vr.x();
        const auto q = vr.y();
        const auto disc = r * r + 4.0 * q;
        if (disc >= 0.0) {
            const auto sqrt_disc = std::sqrt(disc);
            return {{(r + sqrt_disc) / 2.0, 0.0}, {(r - sqrt_disc) / 2.0, 0.0}};
        }
        const auto sqrt_disc = std::sqrt(-disc);
        return {{r / 2.0, sqrt_disc / 2.0}, {r / 2.0, -sqrt_disc / 2.0}};
    }

    /// @brief One Bairstow Newton correction for a quadratic factor (even degree).
    ///
    /// Reads the current factor via `get(idx)`, suppresses all other factors via
    /// the neighbor iterable, and writes the corrected factor via `set(value)`.
    struct even_bairstow_step {
        const std::vector<double>& coeffs;
        std::size_t degree;
        const Options& options;

        template <typename Get, typename Set, typename Neighbors>
        auto operator()(std::size_t idx, Get&& get, Set&& set, Neighbors&& neighbors) const
            -> double {
            const auto vri = get(idx);
            auto local_coeffs = coeffs;  // horner corrupts the array
            auto vA = horner(local_coeffs, degree, vri);
            const auto tol_i = std::max(std::abs(vA.x()), std::abs(vA.y()));
            if (tol_i < options.tol_ind) return 0.0;
            auto vA1 = horner(local_coeffs, degree - 2, vri);
            for (const auto jdx : neighbors) {
                suppress_old(vA, vA1, vri, get(jdx));
            }
            set(vri - delta_scalar(vA, vri, vA1));
            return tol_i;
        }
    };

    /// @brief One Bairstow Newton correction respecting palindromic (autocorr)
    /// symmetry: each neighbor factor contributes both vrj and its reciprocal
    /// image (-vrj.x, 1) / vrj.y; the factor itself also suppresses its own
    /// reciprocal image.
    struct autocorr_bairstow_step {
        const std::vector<double>& coeffs;
        std::size_t degree;
        const Options& options;

        template <typename Get, typename Set, typename Neighbors>
        auto operator()(std::size_t idx, Get&& get, Set&& set, Neighbors&& neighbors) const
            -> double {
            const auto vri = get(idx);
            auto local_coeffs = coeffs;  // horner corrupts the array
            auto vA = horner(local_coeffs, degree, vri);
            const auto tol_i = std::max(std::abs(vA.x()), std::abs(vA.y()));
            if (tol_i < options.tol_ind) return 0.0;
            auto vA1 = horner(local_coeffs, degree - 2, vri);
            for (const auto jdx : neighbors) {
                const auto vrj = get(jdx);
                suppress_old(vA, vA1, vri, vrj);
                const auto vrjn = ginger::Vector2<double>(-vrj.x(), 1.0) / vrj.y();
                suppress_old(vA, vA1, vri, vrjn);
            }
            const auto vrin = ginger::Vector2<double>(-vri.x(), 1.0) / vri.y();
            suppress_old(vA, vA1, vri, vrin);
            set(vri - delta_scalar(vA, vri, vA1));
            return tol_i;
        }
    };

    /// @brief One Aberth-Ehrlich correction for a single root.
    struct aberth_step {
        const std::vector<double>& coeffs;
        const std::vector<double>& coeffs1;  // derivative coefficients of coeffs

        template <typename Get, typename Set, typename Neighbors>
        auto operator()(std::size_t idx, Get&& get, Set&& set, Neighbors&& neighbors) const
            -> double {
            const auto zi = get(idx);
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            for (const auto jdx : neighbors) {
                P1 -= P / (zi - get(jdx));
            }
            set(zi - P / P1);
            return tol_i;
        }
    };

    /// @brief One Aberth-Ehrlich correction for a palindromic (autocorr) root;
    /// each neighbor root contributes both zj and its reciprocal 1/zj.
    struct aberth_autocorr_step {
        const std::vector<double>& coeffs;
        const std::vector<double>& coeffs1;  // derivative coefficients of coeffs

        template <typename Get, typename Set, typename Neighbors>
        auto operator()(std::size_t idx, Get&& get, Set&& set, Neighbors&& neighbors) const
            -> double {
            const auto zi = get(idx);
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            for (const auto jdx : neighbors) {
                P1 -= P / (zi - get(jdx));
                P1 -= P / (zi - 1.0 / get(jdx));
            }
            set(zi - P / P1);
            return tol_i;
        }
    };

    /// @brief Sequential Gauss-Seidel execution: roots are updated in-place in
    /// ascending index order within each iteration.
    struct sequential_policy {
        template <typename State, typename Step>
        static auto run(State& state, const Options& options, Step& step)
            -> std::pair<unsigned int, bool> {
            const auto num_roots = state.size();
            for (auto niter = 0U; niter != options.max_iters; ++niter) {
                auto tolerance = 0.0;
                for (auto idx = std::size_t{0}; idx < num_roots; ++idx) {
                    tolerance = std::max(
                        tolerance,
                        step(
                            idx,
                            [&](std::size_t j) -> typename State::value_type { return state[j]; },
                            [&](typename State::value_type v) { state[idx] = v; },
                            index_range_excluding{num_roots, idx}));
                }
                if (tolerance < options.tolerance) return {niter, true};
            }
            return {options.max_iters, false};
        }
    };

    /// @brief Jacobi multi-threaded execution: each iteration reads a frozen
    /// snapshot of the roots and writes corrected values into disjoint chunks of
    /// the live state, so the iteration is order-independent.  Small problems
    /// (<= 4 roots) fall back to sequential Gauss-Seidel.
    struct jacobi_mt_policy {
        template <typename State, typename Step>
        static auto run(State& state, const Options& options, Step& step)
            -> std::pair<unsigned int, bool> {
            const auto num_roots = state.size();
            if (!ginger::should_parallelize(num_roots)) {
                return sequential_policy::run(state, options, step);
            }
            auto& pool = ginger::get_thread_pool();
            const auto pool_size = pool.size();
            const auto num_threads = std::max(std::size_t{1}, std::min(pool_size, num_roots));
            const auto chunk_size = (num_roots + num_threads - 1) / num_threads;
            for (auto niter = 0U; niter != options.max_iters; ++niter) {
                auto tolerance = 0.0;
                std::vector<std::future<double>> results;
                results.reserve(num_threads);
                auto snapshot = state;
                for (auto t = std::size_t{0}; t < num_threads; ++t) {
                    auto start = t * chunk_size;
                    auto end = std::min(start + chunk_size, num_roots);
                    if (start >= end) break;
                    results.emplace_back(pool.enqueue([&, start, end]() {
                        double max_tol = 0.0;
                        for (auto idx = start; idx < end; ++idx) {
                            max_tol = std::max(
                                max_tol, step(
                                             idx,
                                             [&](std::size_t j) ->
                                             typename State::value_type { return snapshot[j]; },
                                             [&](typename State::value_type v) { state[idx] = v; },
                                             index_range_excluding{num_roots, idx}));
                        }
                        return max_tol;
                    }));
                }
                for (auto&& result : results) {
                    tolerance = std::max(tolerance, result.get());
                }
                if (tolerance < options.tolerance) return {niter, true};
            }
            return {options.max_iters, false};
        }
    };

    /// @brief Atomic decoupled execution: a single `std::atomic` working buffer is
    /// built once and reused across all iterations.  Each thread owns exactly one
    /// slot (single-writer, multi-reader) and iterates independently with no
    /// per-iteration barrier, reading the other slots in round-robin order to
    /// reduce contention.  The returned iteration count is the maximum across
    /// threads and is non-deterministic.
    struct atomic_decoupled_policy {
        template <typename State, typename Step>
        static auto run(State& state, const Options& options, Step& step)
            -> std::pair<unsigned int, bool> {
            using value_type = typename State::value_type;
            const auto num_roots = state.size();
            std::vector<std::atomic<value_type>> buffer(num_roots);
            for (auto idx = std::size_t{0}; idx < num_roots; ++idx) {
                buffer[idx].store(state[idx]);
            }
            const auto use_mt = ginger::should_parallelize(num_roots);
            auto& pool = ginger::get_thread_pool();
            const auto pool_size = pool.size();
            const auto num_threads = use_mt
                                         ? std::max(std::size_t{1}, std::min(pool_size, num_roots))
                                         : std::size_t{1};
            const auto chunk_size
                = use_mt ? (num_roots + num_threads - 1) / num_threads : num_roots;
            fun::Robin<std::size_t> robin(num_roots);
            std::vector<std::future<std::pair<unsigned int, bool>>> results;
            results.reserve(num_threads);
            for (auto t = std::size_t{0}; t < num_threads; ++t) {
                auto start = t * chunk_size;
                auto end = std::min(start + chunk_size, num_roots);
                if (start >= end) break;
                results.emplace_back(pool.enqueue([&, start, end]() {
                    for (auto niter = 0U; niter != options.max_iters; ++niter) {
                        auto max_tol = 0.0;
                        for (auto idx = start; idx < end; ++idx) {
                            max_tol = std::max(
                                max_tol,
                                step(
                                    idx,
                                    [&](std::size_t j) -> value_type { return buffer[j].load(); },
                                    [&](value_type v) { buffer[idx].store(v); },
                                    robin.exclude(idx)));
                        }
                        if (max_tol < options.tolerance) {
                            return std::pair<unsigned int, bool>{niter, true};
                        }
                    }
                    return std::pair<unsigned int, bool>{options.max_iters, false};
                }));
            }
            auto niter = 0U;
            auto converged = true;
            for (auto&& result : results) {
                auto thread_result = result.get();
                niter = std::max(niter, thread_result.first);
                converged = converged && thread_result.second;
            }
            for (auto idx = std::size_t{0}; idx < num_roots; ++idx) {
                state[idx] = buffer[idx].load();
            }
            return {niter, converged};
        }
    };

}  // namespace ginger::detail
