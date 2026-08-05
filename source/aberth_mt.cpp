#include <algorithm>
#include <cmath>    // for acos, cos, sin
#include <complex>  // for complex, operator*, operator+
#include <future>   // for future
#include <ginger/aberth_mt.hpp>
#include <ginger/config.hpp>
#include <ginger/robin.hpp>        // for Robin
#include <ginger/thread_pool.hpp>  // for thread_pool
#include <lds/lds.hpp>
#include <utility>  // for pair
#include <vector>   // for vector, vector<>::reference, __v...

using std::vector;
using Complex = std::complex<double>;

// MT core — uses futures from thread pool with batched scheduling
template <typename F> static auto aberth_mt_core(const vector<double>& coeffs, vector<Complex>& zs,
                                                 const Options& options, ginger::thread_pool& pool,
                                                 F& aberth_job_generator)
    -> std::pair<unsigned int, bool> {
    const auto num_roots = zs.size();
    // For small problems, parallel overhead dominates; run sequentially.
    const auto use_mt = num_roots > 4;
    const auto pool_size = pool.size();
    const auto num_threads
        = use_mt ? std::max(size_t{1}, std::min(pool_size, num_roots)) : size_t{1};
    const auto chunk_size = use_mt ? (num_roots + num_threads - 1) / num_threads : num_roots;

    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        auto zs_snapshot = zs;
        // job reads from zs_snapshot (frozen for the iteration)
        auto aberth_job = aberth_job_generator(coeffs, zs_snapshot);

        vector<std::future<double>> results;
        results.reserve(num_threads);
        for (auto t = size_t{0}; t < num_threads; ++t) {
            auto start = t * chunk_size;
            auto end = std::min(start + chunk_size, num_roots);
            if (start >= end) break;

            results.emplace_back(pool.enqueue([&, start, end]() {
                double max_tol = 0.0;
                for (auto idx = start; idx < end; ++idx) {
                    max_tol = std::max(max_tol, aberth_job(idx));
                }
                return max_tol;
            }));
        }
        for (auto& result : results) {
            tolerance = std::max(tolerance, result.get());
        }
        // copy snapshot updates back to zs
        for (auto idx = size_t{0}; idx < num_roots; ++idx) {
            zs[idx] = zs_snapshot[idx];
        }
        if (tolerance < options.tolerance) {
            return {niter, true};
        }
    }
    return {options.max_iters, false};
}

auto aberth_mt(const vector<double>& coeffs, vector<Complex>& zs,
               const Options& options = Options()) -> std::pair<unsigned int, bool> {
    auto& pool = ginger::get_thread_pool();
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto idx = 0U; idx != degree; ++idx) {
        coeffs1[idx] = static_cast<double>(degree - idx) * coeffs[idx];
    }
    const auto num_zs = zs.size();

    auto aberth_job_generator = [&](const vector<double>&, vector<Complex>& zs_ref) {
        return [&, num_zs](size_t idx) -> double {
            const auto zi = zs_ref[idx];
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            for (auto jdx = 0U; jdx < num_zs; ++jdx) {
                if (jdx == idx) continue;
                P1 -= P / (zi - zs_ref[jdx]);
            }
            zs_ref[idx] -= P / P1;
            return tol_i;
        };
    };

    return aberth_mt_core(coeffs, zs, options, pool, aberth_job_generator);
}

// MT core — uses futures from thread pool
template <typename F>
static auto aberth_autocorr_mt_core(const vector<double>& coeffs, vector<Complex>& zs,
                                    const Options& options, ginger::thread_pool& pool,
                                    F& aberth_job_generator) -> std::pair<unsigned int, bool> {
    const auto num_roots = zs.size();
    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        auto aberth_job = aberth_job_generator(coeffs, zs);
        vector<std::future<double>> results;
        results.reserve(num_roots);
        for (auto idx = 0U; idx != num_roots; ++idx) {
            results.emplace_back(pool.enqueue([&, idx]() { return aberth_job(idx); }));
        }
        for (auto& result : results) {
            tolerance = std::max(tolerance, result.get());
        }
        if (tolerance < options.tolerance) {
            return {niter, true};
        }
    }
    return {options.max_iters, false};
}

auto aberth_autocorr_mt(const vector<double>& coeffs, vector<Complex>& zs,
                        const Options& options = Options()) -> std::pair<unsigned int, bool> {
    auto& pool = ginger::get_thread_pool();
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto idx = 0U; idx != degree; ++idx) {
        coeffs1[idx] = static_cast<double>(degree - idx) * coeffs[idx];
    }
    const auto num_zs = zs.size();

    auto aberth_job_generator = [&](const vector<double>&, vector<Complex>& zs_ref) {
        return [&, num_zs](size_t idx) -> double {
            const auto zi = zs_ref[idx];
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            for (auto jdx = 0U; jdx < num_zs; ++jdx) {
                if (jdx == idx) continue;
                P1 -= P / (zi - zs_ref[jdx]);
                P1 -= P / (zi - 1.0 / zs_ref[jdx]);
            }
            zs_ref[idx] -= P / P1;
            return tol_i;
        };
    };

    return aberth_autocorr_mt_core(coeffs, zs, options, pool, aberth_job_generator);
}
