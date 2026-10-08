/**
 * @file config.hpp
 * @brief Algorithm options (max iterations, tolerance)
 */

#pragma once

#include <algorithm>  // for max
#include <cmath>      // for abs
#include <cstddef>    // for size_t
#include <vector>     // for vector

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

    /**
     * @brief Normalisation factor making the convergence test scale-invariant.
     *
     * At a root, |P(z)| has a floating-point floor proportional to the
     * coefficient magnitudes. A fixed *absolute* tolerance is therefore
     * unreachable for badly scaled polynomials, and the solver burns every
     * iteration up to `max_iters` even though the roots are already accurate.
     * Dividing the coefficients by this factor (roots are invariant under the
     * scaling P -> P/s) makes `Options::tolerance` a relative measure.
     *
     * The factor is floored at 1.0 so the historical absolute test is never
     * made *stricter* for well-scaled polynomials.
     *
     * @param[in] coeffs Polynomial coefficients (highest degree first)
     * @return Scale factor, at least 1.0
     */
    inline auto residual_scale(const std::vector<double>& coeffs) -> double {
        auto scale = 1.0;
        for (const auto coeff : coeffs) {
            scale = std::max(scale, std::abs(coeff));
        }
        return scale;
    }

}  // namespace ginger
