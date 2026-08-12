#include <algorithm>
#include <cmath>              // for abs
#include <cstddef>            // for size_t
#include <future>             // for future
#include <ginger/autocorr_mt.hpp>
#include <ginger/config.hpp>        // for Options
#include <ginger/thread_pool.hpp>   // for get_thread_pool, thread_pool
#include <ginger/vector2.hpp>       // for operator-, Vector2
#include <utility>                  // for pair
#include <vector>                   // for vector

auto pbairstow_autocorr_mt(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                           const Options& options) -> std::pair<unsigned int, bool> {
    auto& pool = ginger::get_thread_pool();
    const auto num_roots = vrs.size();
    const auto degree = coeffs.size() - 1;

    // Small problem: run Gauss-Seidel sequentially (no snapshot overhead)
    if (num_roots <= 4) {
        for (auto niter = 0U; niter != options.max_iters; ++niter) {
            auto tolerance = 0.0;
            for (auto idx = 0U; idx != num_roots; ++idx) {
                const auto& vri = vrs[idx];
                auto local_coeffs = coeffs;  // horner corrupts the array
                auto vA = horner(local_coeffs, degree, vri);
                const auto tol_i = std::max(std::abs(vA.x()), std::abs(vA.y()));
                if (tol_i < options.tol_ind) continue;
                auto vA1 = horner(local_coeffs, degree - 2, vri);
                for (auto jdx = 0U; jdx < num_roots; ++jdx) {
                    if (jdx == idx) continue;
                    const auto& vrj = vrs[jdx];
                    suppress_old(vA, vA1, vri, vrj);
                    const auto vrjn = ginger::Vector2<double>(-vrj.x(), 1.0) / vrj.y();
                    suppress_old(vA, vA1, vri, vrjn);
                }
                const auto vrin = ginger::Vector2<double>(-vri.x(), 1.0) / vri.y();
                suppress_old(vA, vA1, vri, vrin);
                vrs[idx] -= delta_scalar(vA, vri, vA1);
                tolerance = std::max(tolerance, tol_i);
            }
            if (tolerance < options.tolerance) return {niter, true};
        }
        return {options.max_iters, false};
    }

    // Multi-threaded path: batch scheduling + Jacobi snapshot
    const auto pool_size = pool.size();
    const auto num_threads = std::max(size_t{1}, std::min(pool_size, num_roots));
    const auto chunk_size = (num_roots + num_threads - 1) / num_threads;

    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        std::vector<std::future<double>> results;
        results.reserve(num_threads);

        auto vrs_snapshot = vrs;

        for (auto t = size_t{0}; t < num_threads; ++t) {
            auto start = t * chunk_size;
            auto end = std::min(start + chunk_size, num_roots);
            if (start >= end) break;

            results.emplace_back(pool.enqueue(
                [&coeffs, &vrs, &vrs_snapshot, &options, start, end, degree, num_roots]() {
                    double max_tol = 0.0;
                    for (auto idx = start; idx < end; ++idx) {
                        const auto& vri = vrs_snapshot[idx];
                        auto local_coeffs = coeffs;  // horner corrupts the array
                        auto vA = horner(local_coeffs, degree, vri);
                        const auto tol_i = std::max(std::abs(vA.x()), std::abs(vA.y()));
                        if (tol_i < options.tol_ind) continue;
                        auto vA1 = horner(local_coeffs, degree - 2, vri);
                        for (auto jdx = 0U; jdx < num_roots; ++jdx) {
                            if (jdx == idx) continue;
                            const auto& vrj = vrs_snapshot[jdx];
                            suppress_old(vA, vA1, vri, vrj);
                            const auto vrjn = ginger::Vector2<double>(-vrj.x(), 1.0) / vrj.y();
                            suppress_old(vA, vA1, vri, vrjn);
                        }
                        const auto vrin = ginger::Vector2<double>(-vri.x(), 1.0) / vri.y();
                        suppress_old(vA, vA1, vri, vrin);

                        vrs[idx] -= delta_scalar(vA, vri, vA1);
                        max_tol = std::max(max_tol, tol_i);
                    }
                    return max_tol;
                }));
        }
        for (auto&& result : results) {
            auto&& res = result.get();
            tolerance = std::max(tolerance, res);
        }
        if (tolerance < options.tolerance) {
            return {niter, true};
        }
    }
    return {options.max_iters, false};
}
