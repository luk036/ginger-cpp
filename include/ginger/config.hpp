/**
 * @file config.hpp
 * @brief Algorithm options (max iterations, tolerance)
 */

#pragma once

#include <cstddef>  // for size_t

/**
 * @brief Options for convergence-based algorithms
 *
 * Stores maximum iteration count, tolerance for global convergence,
 * and per-root tolerance for convergence checks
 * used by Bairstow and Aberth root-finding methods.
 */
namespace ginger {

    class Options {
      public:
        unsigned int max_iters = 2000U;
        double tolerance = 1e-12;
        double tol_ind = 1e-15;
    };

    /// @brief Number of roots above which the multi-threaded policies are used by default.
    inline constexpr std::size_t PARALLEL_THRESHOLD = 4;

    /// @brief Whether `num_roots` should use the multi-threaded execution policy.
    inline auto should_parallelize(std::size_t num_roots) -> bool {
        return num_roots > PARALLEL_THRESHOLD;
    }

}  // namespace ginger
