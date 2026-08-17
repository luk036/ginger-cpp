#include <algorithm>
#include <atomic>                         // for atomic
#include <cmath>                          // for abs
#include <cstddef>                        // for size_t
#include <future>                         // for future
#include <ginger/config.hpp>              // for Options
#include <ginger/robin.hpp>               // for Robin
#include <ginger/rootfinding_atomic.hpp>  // for Vec2, delta_scalar, horner, pbairstow_even_atomic
#include <ginger/thread_pool.hpp>         // for get_thread_pool, thread_pool
#include <utility>                        // for pair
#include <vector>                         // for vector

auto pbairstow_even_atomic(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                           const Options& options) -> std::pair<unsigned int, bool> {
    auto& pool = ginger::get_thread_pool();
    const auto num_roots = vrs.size();
    const auto degree = coeffs.size() - 1;

    // Atomic working buffer built ONCE (whole-pair atomicity: correctness over lock-free).
    std::vector<std::atomic<Vec2>> vrs_atomic(num_roots);
    for (auto idx = size_t{0}; idx < num_roots; ++idx) {
        vrs_atomic[idx].store(vrs[idx]);
    }

    // For small problems, parallel overhead dominates; run single-threaded.
    const auto use_mt = num_roots > 4;
    const auto pool_size = pool.size();
    const auto num_threads
        = use_mt ? std::max(size_t{1}, std::min(pool_size, num_roots)) : size_t{1};
    const auto chunk_size = use_mt ? (num_roots + num_threads - 1) / num_threads : num_roots;
    fun::Robin<size_t> robin(num_roots);

    // Decoupled: each thread runs its own iteration loop INDEPENDENTLY (no per-iteration
    // barrier). A thread exits as soon as its own chunk converges or max_iters is exceeded.
    std::vector<std::future<std::pair<unsigned int, bool>>> results;
    results.reserve(num_threads);
    for (auto t = size_t{0}; t < num_threads; ++t) {
        auto start = t * chunk_size;
        auto end = std::min(start + chunk_size, num_roots);
        if (start >= end) break;

        results.emplace_back(pool.enqueue([&, start, end]() {
            for (auto niter = 0U; niter != options.max_iters; ++niter) {
                auto max_tol = 0.0;
                for (auto idx = start; idx < end; ++idx) {
                    const auto vri = vrs_atomic[idx].load();
                    auto local_coeffs = coeffs;  // horner corrupts the array
                    auto vA = horner(local_coeffs, degree, vri);
                    const auto tol_i = std::max(std::abs(vA.x()), std::abs(vA.y()));
                    if (tol_i < options.tol_ind) continue;
                    auto vA1 = horner(local_coeffs, degree - 2, vri);
                    // Round-robin suppression order: each thread reads the other
                    // slots in a different rotation, reducing concurrent access
                    // to the same atomic variable.
                    for (auto jdx : robin.exclude(idx)) {
                        const auto vrj = vrs_atomic[jdx].load();
                        suppress_old(vA, vA1, vri, vrj);
                    }
                    vrs_atomic[idx].store(vri - delta_scalar(vA, vri, vA1));
                    max_tol = std::max(max_tol, tol_i);
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

    // Publish the final atomic buffer back to the caller.
    for (auto idx = size_t{0}; idx < num_roots; ++idx) {
        vrs[idx] = vrs_atomic[idx].load();
    }
    return {niter, converged};
}
