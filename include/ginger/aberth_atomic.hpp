/**
 * @file aberth_atomic.hpp
 * @brief Aberth-Ehrlich method for polynomial root-finding (atomic)
 */

#pragma once

#include "aberth.hpp"

/**
 * @brief Atomic Aberth-Ehrlich method
 *
 * The `aberth_atomic` function is an implementation of the Aberth-Ehrlich method for finding the
 * roots of a polynomial using an atomic working buffer.
 *
 * Atomic variant of the Aberth-Ehrlich method:
 * @f[
 *     x_k^{(i+1)} = x_k^{(i)} - \frac{P(x_k)}{P'(x_k)}\Bigg/ \left(1 -
 * \frac{P(x_k)}{P'(x_k)}\sum_{j \ne k}\frac{1}{x_k - x_j}\right)
 * @f]
 * The atomic working buffer is built once (no per-iteration snapshots). Each thread owns exactly
 * one slot (single-writer, multi-reader): jobs `load()` the other slots and `store()` only their
 * own slot, so the iteration becomes asynchronous and in-place (Gauss-Seidel-like). Threads run
 * independently with no per-iteration synchronization: each thread exits when its own roots
 * converge or the maximum number of iterations is exceeded. The whole
 * (re, im) pair is read and written atomically via `std::atomic<std::complex<double>>`
 * (16 bytes, not lock-free on MSVC but correct).
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of a
 * polynomial. Each element of the vector corresponds to a term in the polynomial, starting from the
 * highest degree term and ending with the constant term. For example, if the polynomial is
 * `3x^2 + 2x + 1`, the coefficients vector would be `{3, 2, 1}`.
 * @param[in,out] zs `zs` is a vector of complex numbers representing the initial guesses for the
 * roots of the polynomial. The function will update these values iteratively to converge to the
 * actual roots.
 * @param[in] options The `options` parameter is an object of type `Options` that contains the
 * maximum number of iterations (`max_iters`) and the tolerance (`tolerance`). These options control
 * the convergence criteria for the Aberth-Ehrlich method.
 *
 * @return The `aberth_atomic` function returns a `std::pair<unsigned int, bool>`. The first element
 * of the pair represents the number of iterations performed, and the second element represents
 * whether the method converged to a solution within the specified tolerance. Iteration counts are
 * non-deterministic (each thread may take a different number of iterations; the returned count is
 * the maximum across threads), so tests should assert convergence, not a fixed count.
 */
extern auto aberth_atomic(const std::vector<double>& coeffs, std::vector<std::complex<double>>& zs,
                          const Options& options) -> std::pair<unsigned int, bool>;

/**
 * @brief Atomic Aberth-Ehrlich method (specifically for auto-correlation functions)
 *
 * The `aberth_autocorr_atomic` function is an implementation of the Aberth-Ehrlich method for
 * finding the roots of a palindromic (auto-correlation) polynomial using an atomic working buffer.
 *
 * The atomic working buffer is built once (no per-iteration snapshots). Each thread owns exactly
 * one slot (single-writer, multi-reader): jobs `load()` the other slots and `store()` only their
 * own slot, so the iteration becomes asynchronous and in-place (Gauss-Seidel-like). Threads run
 * independently with no per-iteration synchronization: each thread exits when its own roots
 * converge or the maximum number of iterations is exceeded. The whole
 * (re, im) pair is read and written atomically via `std::atomic<std::complex<double>>`
 * (16 bytes, not lock-free on MSVC but correct).
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of a
 * polynomial. Each element of the vector corresponds to a term in the polynomial, starting from the
 * highest degree term and ending with the constant term. For example, if the polynomial is
 * `3x^2 + 2x + 1`, the coefficients vector would be `{3, 2, 1}`.
 * @param[in,out] zs `zs` is a vector of complex numbers representing the initial guesses for the
 * roots of the polynomial. The function will update these values iteratively to converge to the
 * actual roots.
 * @param[in] options The `options` parameter is an object of type `Options` that contains the
 * maximum number of iterations (`max_iters`) and the tolerance (`tolerance`). These options control
 * the convergence criteria for the Aberth-Ehrlich method.
 *
 * @return The `aberth_autocorr_atomic` function returns a `std::pair<unsigned int, bool>`. The
 * first element of the pair represents the number of iterations performed, and the second element
 * represents whether the method converged to a solution within the specified tolerance. Iteration
 * counts are non-deterministic (each thread may take a different number of iterations; the returned
 * count is the maximum across threads), so tests should assert convergence, not a fixed count.
 */
extern auto aberth_autocorr_atomic(const std::vector<double>& coeffs,
                                   std::vector<std::complex<double>>& zs, const Options& options)
    -> std::pair<unsigned int, bool>;
