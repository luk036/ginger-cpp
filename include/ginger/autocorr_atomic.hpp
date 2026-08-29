/**
 * @file autocorr_atomic.hpp
 * @brief Auto-correlation polynomial root-finding (palindromic polynomials) (atomic)
 */

#pragma once

#include "autocorr.hpp"

/**
 * @brief Atomic multi-threading Bairstow's method (specific for auto-correlation function)
 *
 * The function `pbairstow_autocorr_atomic` implements Bairstow's method for finding the roots of a
 * palindromic (auto-correlation) polynomial using a single atomic working buffer.
 *
 * Unlike the snapshot-based Jacobi variant (`pbairstow_autocorr_mt`), a single atomic working
 * buffer is built once and reused across all iterations: each thread owns exactly one slot
 * (single-writer, multi-reader) and reads the latest values of the other slots while iterating.
 * The whole @f$ (r_i, q_i) @f$ pair is read/written atomically via `std::atomic<Vec2>` (16 bytes,
 * not lock-free on MSVC but correct). Iterations are asynchronous/in-place (Gauss-Seidel-like):
 * each thread runs its own iteration loop independently with no per-iteration synchronization
 * and exits when its own factors converge or the maximum number of iterations is exceeded. The
 * iteration count is therefore NON-DETERMINISTIC (the returned count is the maximum across
 * threads) — convergence, not a fixed iteration count, must be asserted.
 *
 * Each thread finds a quadratic factor @f$ x^2 - r_i x - q_i @f$ that respects the palindromic
 * symmetry, where the roots appear in reciprocal pairs:
 * @f[
 *     \begin{bmatrix} \Delta r_i \\ \Delta q_i \end{bmatrix} = -J_i^{-1} \begin{bmatrix} P_i \\ Q_i
 * \end{bmatrix}
 * @f]
 *
 * @param[in] coeffs polynomial
 * @param[in,out] vrs vector of iterates
 * @param[in] options maximum iterations and tolorance
 * @return std::pair<unsigned int, bool>
 */
extern auto pbairstow_autocorr_atomic(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                                      const ginger::Options& options)
    -> std::pair<unsigned int, bool>;
