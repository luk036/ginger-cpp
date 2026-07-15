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
 *
 * @f$ \phi_2(n) = \sum_{k=0}^{\infty} a_k(n) \, 2^{-k-1} @f$
 * where @f$ a_k(n) @f$ are the binary digits of @f$ n @f$.
 *
 * @param[in] index Sequence index
 * @return double The van der Corput value
 */
extern double vdc2_table(unsigned long index);

/**
 * @brief Cosine of pi times van der Corput value
 *
 * @f$ \cos(\pi \cdot \phi_2(\text{index})) @f$
 *
 * @param[in] index Sequence index
 * @return double cos(pi * vdc2_table(index))
 */
extern double cos_pi_vdc2(unsigned long index);

/**
 * @brief Circle sequence X-coordinate for an index
 *
 * @f$ x = \cos(2\pi \cdot \phi_2(\text{index})) @f$
 *
 * @param[in] index Sequence index
 * @return double X-coordinate on the unit circle
 */
extern double circle2_table_x(unsigned long index);

/**
 * @brief Circle sequence Y-coordinate for an index
 *
 * @f$ y = \sin(2\pi \cdot \phi_2(\text{index})) @f$
 *
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
 * @f$ z_k = R \cdot e^{2\pi i \cdot \phi_2(k)}, \quad k = 0,\dots,n-1 @f$
 * where @f$ R @f$ is estimated from the polynomial coefficients.
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
 * @dot
 *   digraph aberth_iter {
 *     rankdir=LR; bgcolor="transparent";
 *     node [shape=box, style=filled, fillcolor="#d4e6f1"];
 *     init [label="Initial\nguesses z_k", fillcolor="#a9cce3"];
 *     eval [label="Evaluate\nP(z_k), P'(z_k)"];
 *     update [label="Update\nz_k_new"];
 *     check [label="Converged?", shape=diamond, fillcolor="#f9e79f"];
 *     done [label="Roots found!", fillcolor="#7fb3d8"];
 *     init -> eval -> update -> check;
 *     check -> eval [label="No", style=dashed, color="#e74c3c"];
 *     check -> done [label="Yes", color="#27ae60"];
 *   }
 * @enddot
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
 * The `initial_aberth_autocorr` function calculates the initial values for the Aberth-Ehrlich
 * method for finding the roots of a palindromic (auto-correlation) polynomial.
 *
 * @f$ z_k = e^{2\pi i \cdot \phi_2(k)}, \quad k = 0,\dots,\lfloor n/2 \rfloor @f$
 * where the roots are distributed on the unit circle.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector of doubles that represents the coefficients
 * of a polynomial.
 *
 * @return The function `initial_aberth_autocorr` returns a vector of Complex numbers.
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
 * @f[
 *     p_0 = \arg\min_{z \in S} |z|, \qquad
 *     p_k = \arg\max_{z \in S_k} \min_{j < k} |z - p_j|
 * @f]
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
 * @f[
 *     P(x) = \prod_{k=1}^{n} (x - r_k)(x - r_k^{-1}) = x^{2n} + a_{2n-1}x^{2n-1} + \cdots + a_0
 * @f]
 *
 * @param[in] zs Roots found by aberth_autocorr or aberth_autocorr_mt
 * @return std::vector<double> Monic polynomial coefficients (highest degree first)
 */
extern auto poly_from_autocorr_roots(const std::vector<std::complex<double>>& zs)
    -> std::vector<double>;
