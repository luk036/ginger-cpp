/**
 * @file solve.hpp
 * @brief Facade entry points that auto-select the execution policy
 */

#pragma once

#include <complex>  // for complex
#include <utility>  // for pair
#include <vector>   // for vector

#include "aberth.hpp"              // for aberth, aberth_autocorr
#include "aberth_atomic.hpp"       // for aberth_atomic, aberth_autocorr_atomic
#include "aberth_mt.hpp"           // for aberth_mt, aberth_autocorr_mt
#include "autocorr.hpp"            // for pbairstow_autocorr_st
#include "autocorr_atomic.hpp"     // for pbairstow_autocorr_atomic
#include "autocorr_mt.hpp"         // for pbairstow_autocorr_mt
#include "config.hpp"              // for ginger::Options, should_parallelize
#include "rootfinding.hpp"         // for pbairstow_even_st
#include "rootfinding_atomic.hpp"  // for pbairstow_even_atomic
#include "rootfinding_mt.hpp"      // for pbairstow_even_mt

namespace ginger {

    /// @brief Execution policy selection for the solver facades.
    enum class solve_mode { automatic, sequential, multi_threaded, atomic };

    /// @brief Facade for the parallel Bairstow method (even degree).
    ///
    /// Selects the execution policy: `automatic` (default) dispatches to the
    /// multi-threaded variant when `should_parallelize(vrs.size())` holds,
    /// otherwise to the single-threaded variant.
    ///
    /// @param[in] coeffs Polynomial coefficients (highest degree first)
    /// @param[in,out] vrs Quadratic factor iterates
    /// @param[in] options Convergence options
    /// @param[in] mode Execution policy selector
    /// @return std::pair<unsigned int, bool> (iterations, converged)
    inline auto solve_pbairstow_even(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                                     const Options& options,
                                     solve_mode mode = solve_mode::automatic)
        -> std::pair<unsigned int, bool> {
        switch (mode) {
            case solve_mode::sequential:
                return pbairstow_even_st(coeffs, vrs, options);
            case solve_mode::multi_threaded:
                return pbairstow_even_mt(coeffs, vrs, options);
            case solve_mode::atomic:
                return pbairstow_even_atomic(coeffs, vrs, options);
            case solve_mode::automatic:
            default:
                return should_parallelize(vrs.size()) ? pbairstow_even_mt(coeffs, vrs, options)
                                                      : pbairstow_even_st(coeffs, vrs, options);
        }
    }

    /// @brief Facade for the parallel Bairstow method (auto-correlation).
    /// @param[in] coeffs Polynomial coefficients (highest degree first)
    /// @param[in,out] vrs Quadratic factor iterates
    /// @param[in] options Convergence options
    /// @param[in] mode Execution policy selector
    /// @return std::pair<unsigned int, bool> (iterations, converged)
    inline auto solve_pbairstow_autocorr(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                                         const Options& options,
                                         solve_mode mode = solve_mode::automatic)
        -> std::pair<unsigned int, bool> {
        switch (mode) {
            case solve_mode::sequential:
                return pbairstow_autocorr_st(coeffs, vrs, options);
            case solve_mode::multi_threaded:
                return pbairstow_autocorr_mt(coeffs, vrs, options);
            case solve_mode::atomic:
                return pbairstow_autocorr_atomic(coeffs, vrs, options);
            case solve_mode::automatic:
            default:
                return should_parallelize(vrs.size()) ? pbairstow_autocorr_mt(coeffs, vrs, options)
                                                      : pbairstow_autocorr_st(coeffs, vrs, options);
        }
    }

    /// @brief Facade for the Aberth-Ehrlich method.
    /// @param[in] coeffs Polynomial coefficients (highest degree first)
    /// @param[in,out] zs Root iterates
    /// @param[in] options Convergence options
    /// @param[in] mode Execution policy selector
    /// @return std::pair<unsigned int, bool> (iterations, converged)
    inline auto solve_aberth(const std::vector<double>& coeffs,
                             std::vector<std::complex<double>>& zs, const Options& options,
                             solve_mode mode = solve_mode::automatic)
        -> std::pair<unsigned int, bool> {
        switch (mode) {
            case solve_mode::sequential:
                return aberth(coeffs, zs, options);
            case solve_mode::multi_threaded:
                return aberth_mt(coeffs, zs, options);
            case solve_mode::atomic:
                return aberth_atomic(coeffs, zs, options);
            case solve_mode::automatic:
            default:
                return should_parallelize(zs.size()) ? aberth_mt(coeffs, zs, options)
                                                     : aberth(coeffs, zs, options);
        }
    }

    /// @brief Facade for the Aberth-Ehrlich method (auto-correlation).
    /// @param[in] coeffs Polynomial coefficients (highest degree first)
    /// @param[in,out] zs Root iterates
    /// @param[in] options Convergence options
    /// @param[in] mode Execution policy selector
    /// @return std::pair<unsigned int, bool> (iterations, converged)
    inline auto solve_aberth_autocorr(const std::vector<double>& coeffs,
                                      std::vector<std::complex<double>>& zs, const Options& options,
                                      solve_mode mode = solve_mode::automatic)
        -> std::pair<unsigned int, bool> {
        switch (mode) {
            case solve_mode::sequential:
                return aberth_autocorr(coeffs, zs, options);
            case solve_mode::multi_threaded:
                return aberth_autocorr_mt(coeffs, zs, options);
            case solve_mode::atomic:
                return aberth_autocorr_atomic(coeffs, zs, options);
            case solve_mode::automatic:
            default:
                return should_parallelize(zs.size()) ? aberth_autocorr_mt(coeffs, zs, options)
                                                     : aberth_autocorr(coeffs, zs, options);
        }
    }

}  // namespace ginger
