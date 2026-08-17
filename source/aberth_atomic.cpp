#include <algorithm>
#include <atomic>   // for atomic
#include <cmath>    // for abs
#include <complex>  // for complex, operator*, operator+
#include <future>   // for future
#include <ginger/aberth_atomic.hpp>
#include <ginger/config.hpp>
#include <ginger/robin.hpp>        // for Robin
#include <ginger/thread_pool.hpp>  // for thread_pool
#include <lds/lds.hpp>
#include <utility>  // for pair
#include <vector>   // for vector

using std::vector;
using Complex = std::complex<double>;
using AtomicComplex = std::atomic<Complex>;

// Atomic core: one atomic working buffer built ONCE. Each job reads all slots
// via .load() and writes only its own slot via .store() (single-writer,
// multi-reader) -> asynchronous in-place updates, no per-iteration snapshots.
template <typename F> static auto aberth_atomic_core(const vector<double>& coeffs,
                                                     vector<Complex>& zs, const Options& options,
                                                     ginger::thread_pool& pool, F& job_generator)
    -> std::pair<unsigned int, bool> {
    const auto num_roots = zs.size();
    // Build the atomic working buffer ONCE (whole-pair atomics: correctness over speed).
    vector<AtomicComplex> zs_atomic(num_roots);
    for (auto idx = size_t{0}; idx < num_roots; ++idx) {
        zs_atomic[idx].store(zs[idx]);
    }
    auto job = job_generator(coeffs, zs_atomic);

    // For small problems, parallel overhead dominates; run sequentially.
    const auto use_mt = num_roots > 4;
    const auto pool_size = pool.size();
    const auto num_threads
        = use_mt ? std::max(size_t{1}, std::min(pool_size, num_roots)) : size_t{1};
    const auto chunk_size = use_mt ? (num_roots + num_threads - 1) / num_threads : num_roots;

    // Decoupled: each thread runs its own iteration loop INDEPENDENTLY (no per-iteration
    // barrier). A thread exits as soon as its own chunk converges or max_iters is exceeded.
    vector<std::future<std::pair<unsigned int, bool>>> results;
    results.reserve(num_threads);
    for (auto t = size_t{0}; t < num_threads; ++t) {
        auto start = t * chunk_size;
        auto end = std::min(start + chunk_size, num_roots);
        if (start >= end) break;
        results.emplace_back(pool.enqueue([&, start, end]() {
            for (auto niter = 0U; niter != options.max_iters; ++niter) {
                auto max_tol = 0.0;
                for (auto idx = start; idx < end; ++idx) {
                    max_tol = std::max(max_tol, job(idx));
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
    for (auto& result : results) {
        auto thread_result = result.get();
        niter = std::max(niter, thread_result.first);
        converged = converged && thread_result.second;
    }

    // Publish the final atomic buffer back to the caller.
    for (auto idx = size_t{0}; idx < num_roots; ++idx) {
        zs[idx] = zs_atomic[idx].load();
    }
    return {niter, converged};
}

auto aberth_atomic(const vector<double>& coeffs, vector<Complex>& zs,
                   const Options& options = Options()) -> std::pair<unsigned int, bool> {
    auto& pool = ginger::get_thread_pool();
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto idx = 0U; idx != degree; ++idx) {
        coeffs1[idx] = static_cast<double>(degree - idx) * coeffs[idx];
    }
    const auto num_zs = zs.size();
    fun::Robin<size_t> robin(num_zs);

    auto aberth_job_generator = [&](const vector<double>&, vector<AtomicComplex>& zs_ref) {
        return [&, num_zs](size_t idx) -> double {
            const auto zi = zs_ref[idx].load();
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            // Round-robin suppression order: each thread reads the other slots
            // in a different rotation, reducing concurrent access to the same slot.
            for (auto jdx : robin.exclude(idx)) {
                P1 -= P / (zi - zs_ref[jdx].load());
            }
            zs_ref[idx].store(zi - P / P1);
            return tol_i;
        };
    };

    return aberth_atomic_core(coeffs, zs, options, pool, aberth_job_generator);
}

auto aberth_autocorr_atomic(const vector<double>& coeffs, vector<Complex>& zs,
                            const Options& options = Options()) -> std::pair<unsigned int, bool> {
    auto& pool = ginger::get_thread_pool();
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto idx = 0U; idx != degree; ++idx) {
        coeffs1[idx] = static_cast<double>(degree - idx) * coeffs[idx];
    }
    const auto num_zs = zs.size();
    fun::Robin<size_t> robin(num_zs);

    auto aberth_job_generator = [&](const vector<double>&, vector<AtomicComplex>& zs_ref) {
        return [&, num_zs](size_t idx) -> double {
            const auto zi = zs_ref[idx].load();
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            // Round-robin suppression order: each thread reads the other slots
            // in a different rotation, reducing concurrent access to the same slot.
            for (auto jdx : robin.exclude(idx)) {
                P1 -= P / (zi - zs_ref[jdx].load());
                P1 -= P / (zi - 1.0 / zs_ref[jdx].load());
            }
            zs_ref[idx].store(zi - P / P1);
            return tol_i;
        };
    };

    return aberth_atomic_core(coeffs, zs, options, pool, aberth_job_generator);
}
