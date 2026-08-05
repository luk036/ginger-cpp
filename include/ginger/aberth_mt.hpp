/**
 * @file aberth_mt.hpp
 * @brief Aberth-Ehrlich method for polynomial root-finding (MT)
 */

#pragma once

#include "aberth.hpp"

/**
 * @brief Multi-threading Aberth-Ehrlich method
 *
 * The `aberth_mt` function is a multi-threaded implementation of the Aberth-Ehrlich method for
 * finding the roots of a polynomial.
 *
 * Multi-threaded variant of the Aberth-Ehrlich method:
 * @f[
 *     x_k^{(i+1)} = x_k^{(i)} - \frac{P(x_k)}{P'(x_k)}\Bigg/ \left(1 -
 * \frac{P(x_k)}{P'(x_k)}\sum_{j \ne k}\frac{1}{x_k - x_j}\right)
 * @f]
 * Each root is updated in parallel using separate threads.
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
 * @return The `aberth_mt` function returns a `std::pair<unsigned int, bool>`. The first element of
 * the pair represents the number of iterations performed, and the second element represents whether
 * the method converged to a solution within the specified tolerance.
 */
extern auto aberth_mt(const std::vector<double>& coeffs, std::vector<std::complex<double>>& zs,
                      const Options& options) -> std::pair<unsigned int, bool>;

/**
 * @brief Multi-threading Aberth-Ehrlich method (specifically for auto-correlation functions)
 *
 * The `aberth_autocorr_mt` function is a multi-threaded implementation of the Aberth-Ehrlich method
 * for finding the roots of a palindromic (auto-correlation) polynomial.
 *
 * Aberth's method is a method for finding the roots of a polynomial that is
 * robust but requires complex arithmetic even if the polynomial is real. This
 * is because it starts with complex initial approximations.
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
 * @return The `aberth_autocorr_mt` function returns a `std::pair<unsigned int, bool>`. The first
 * element of the pair represents the number of iterations performed, and the second element
 * represents whether the method converged to a solution within the specified tolerance.
 */
extern auto aberth_autocorr_mt(const std::vector<double>& coeffs,
                               std::vector<std::complex<double>>& zs, const Options& options)
    -> std::pair<unsigned int, bool>;
