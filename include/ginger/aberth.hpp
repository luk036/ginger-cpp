/**
 * @file aberth.hpp
 * @brief Aberth-Ehrlich method for polynomial root-finding
 */

#pragma once

// import numpy as np
#include <complex>
#include <utility>
#include <vector>

class Options;

/**
 * @brief van der Corput sequence value for a given index
 * @param[in] index Sequence index
 * @return double The van der Corput value
 */
extern double vdc2_table(unsigned long index);

/**
 * @brief Cosine of pi times van der Corput value
 * @param[in] index Sequence index
 * @return double cos(pi * vdc2_table(index))
 */
extern double cos_pi_vdc2(unsigned long index);

/**
 * @brief Circle sequence X-coordinate for an index
 * @param[in] index Sequence index
 * @return double X-coordinate on the unit circle
 */
extern double circle2_table_x(unsigned long index);

/**
 * @brief Circle sequence Y-coordinate for an index
 * @param[in] index Sequence index
 * @return double Y-coordinate on the unit circle
 */
extern double circle2_table_y(unsigned long index);

/**
 * @brief Initial guess for the Aberth-Ehrlich method
 *
 * The `initial_aberth` function calculates the initial values for the Aberth-Ehrlich method for
 * finding the roots of a polynomial.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector of doubles that represents the coefficients
 * of a polynomial.
 *
 * @return The function `initial_aberth` returns a vector of Complex numbers.
 */
extern auto initial_aberth(const std::vector<double>& coeffs) -> std::vector<std::complex<double>>;

/**
 * @brief Single-threading Aberth-Ehrlich method
 *
 * The `aberth` function is an implementation of the Aberth-Ehrlich method for finding
 * the roots of a polynomial.
 *
 * Aberth's method iteratively improves root estimates for a polynomial \f$P(x)\f$:
 * @f[
 *     x_k^{(i+1)} = x_k^{(i)} - \frac{P(x_k)}{P'(x_k)}\Bigg/ \left(1 -
 * \frac{P(x_k)}{P'(x_k)}\sum_{j \ne k}\frac{1}{x_k - x_j}\right)
 * @f]
 * where the sum is over all other root approximations. The method is
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
 * @return The `aberth` function returns a `std::pair<unsigned int, bool>`. The first element of the
 * pair represents the number of iterations performed, and the second element represents whether the
 * method converged to a solution within the specified tolerance.
 */
extern auto aberth(const std::vector<double>& coeffs, std::vector<std::complex<double>>& zs,
                   const Options& options) -> std::pair<unsigned int, bool>;

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
 * @brief Initial guess for the Aberth-Ehrlich method (specifically for auto-correlation functions)
 *
 * The `initial_aberth` function calculates the initial values for the Aberth-Ehrlich method for
 * finding the roots of a polynomial.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector of doubles that represents the coefficients
 * of a polynomial.
 *
 * @return The function `initial_aberth` returns a vector of Complex numbers.
 */
extern auto initial_aberth_autocorr(const std::vector<double>& coeffs)
    -> std::vector<std::complex<double>>;

/**
 * @brief Single-threading Aberth-Ehrlich method (specifically for auto-correlation functions)
 *
 * The `aberth_autocorr` function is an implementation of the Aberth-Ehrlich method for finding
 * the roots of a palindromic (auto-correlation) polynomial.
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
 * @return The `aberth_autocorr` function returns a `std::pair<unsigned int, bool>`. The first
 * element of the pair represents the number of iterations performed, and the second element
 * represents whether the method converged to a solution within the specified tolerance.
 */
extern auto aberth_autocorr(const std::vector<double>& coeffs,
                            std::vector<std::complex<double>>& zs, const Options& options)
    -> std::pair<unsigned int, bool>;

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

/**
 * @brief Reconstruct a monic polynomial from its complex roots
 *
 * Builds the monic polynomial whose roots are the given complex numbers:
 * @f[
 *     P(x) = \prod_{k=1}^{n} (x - z_k) = x^n + a_{n-1}x^{n-1} + \cdots + a_0
 * @f]
 * Uses Leja ordering for numerical stability.
 *
 * @param[in] zs Vector of complex roots
 * @return std::vector<double> Monic polynomial coefficients (highest degree first)
 */
extern auto poly_from_roots(const std::vector<std::complex<double>>& zs) -> std::vector<double>;

/**
 * @brief Leja ordering of complex points
 *
 * Reorders complex points using the greedy Leja algorithm: starts with the
 * smallest-magnitude point, then iteratively selects the remaining point that
 * maximizes the minimum Euclidean distance to all already-selected points.
 * This ordering reduces numerical error when reconstructing polynomials from roots.
 *
 * @param[in] points Input vector of complex numbers
 * @return std::vector<std::complex<double>> Reordered points in Leja sequence
 */
extern auto leja_order(const std::vector<std::complex<double>>& points)
    -> std::vector<std::complex<double>>;

/**
 * @brief Reconstruct a monic polynomial from its autocorrelation roots
 *
 * Auto-correlation (palindromic) polynomials have roots in reciprocal pairs.
 * The aberth_autocorr functions find the degree/2 "independent" roots.
 * This function adds the reciprocal of each root (1/z) to get the full set
 * of degree roots, then reconstructs with Leja ordering.
 *
 * @param[in] zs Roots found by aberth_autocorr or aberth_autocorr_mt
 * @return std::vector<double> Monic polynomial coefficients (highest degree first)
 */
extern auto poly_from_autocorr_roots(const std::vector<std::complex<double>>& zs)
    -> std::vector<double>;
