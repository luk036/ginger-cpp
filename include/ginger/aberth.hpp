#pragma once

// import numpy as np
#include <complex>
#include <utility>
#include <vector>

class Options;

extern double vdc2_table(unsigned long index);
extern double cos_pi_vdc2(unsigned long index);
extern double circle2_table_x(unsigned long index);
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
 * The `aberth` function is a multi-threaded implementation of the Aberth-Ehrlich method for finding
 * the roots of a polynomial.
 *
 * Aberth's method is a method for finding the roots of a polynomial that is
 * robust but requires complex arithmetic even if the polynomial is real. This
 * is because it starts with complex initial approximations.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of a
 * polynomial. Each element of the vector corresponds to a term in the polynomial, starting from the
 * highest degree term and ending with the constant term. For example, if the polynomial is `3x^2 +
 * 2x +
 * @param[in] zs `zs` is a vector of complex numbers representing the initial guesses for the roots
 * of the polynomial. The function will update these values iteratively to converge to the actual
 * roots.
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
 * The `aberth` function is a multi-threaded implementation of the Aberth-Ehrlich method for finding
 * the roots of a polynomial.
 *
 * Aberth's method is a method for finding the roots of a polynomial that is
 * robust but requires complex arithmetic even if the polynomial is real. This
 * is because it starts with complex initial approximations.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of a
 * polynomial. Each element of the vector corresponds to a term in the polynomial, starting from the
 * highest degree term and ending with the constant term. For example, if the polynomial is `3x^2 +
 * 2x +
 * @param[in] zs `zs` is a vector of complex numbers representing the initial guesses for the roots
 * of the polynomial. The function will update these values iteratively to converge to the actual
 * roots.
 * @param[in] options The `options` parameter is an object of type `Options` that contains the
 * maximum number of iterations (`max_iters`) and the tolerance (`tolerance`). These options control
 * the convergence criteria for the Aberth-Ehrlich method.
 *
 * @return The `aberth` function returns a `std::pair<unsigned int, bool>`. The first element of the
 * pair represents the number of iterations performed, and the second element represents whether the
 * method converged to a solution within the specified tolerance.
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
 * The `aberth` function is a multi-threaded implementation of the Aberth-Ehrlich method for finding
 * the roots of a polynomial.
 *
 * Aberth's method is a method for finding the roots of a polynomial that is
 * robust but requires complex arithmetic even if the polynomial is real. This
 * is because it starts with complex initial approximations.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of a
 * polynomial. Each element of the vector corresponds to a term in the polynomial, starting from the
 * highest degree term and ending with the constant term. For example, if the polynomial is `3x^2 +
 * 2x +
 * @param[in] zs `zs` is a vector of complex numbers representing the initial guesses for the roots
 * of the polynomial. The function will update these values iteratively to converge to the actual
 * roots.
 * @param[in] options The `options` parameter is an object of type `Options` that contains the
 * maximum number of iterations (`max_iters`) and the tolerance (`tolerance`). These options control
 * the convergence criteria for the Aberth-Ehrlich method.
 *
 * @return The `aberth` function returns a `std::pair<unsigned int, bool>`. The first element of the
 * pair represents the number of iterations performed, and the second element represents whether the
 * method converged to a solution within the specified tolerance.
 */
extern auto aberth_autocorr(const std::vector<double>& coeffs,
                            std::vector<std::complex<double>>& zs, const Options& options)
    -> std::pair<unsigned int, bool>;

/**
 * @brief Multi-threading Aberth-Ehrlich method (specifically for auto-correlation functions)
 *
 * The `aberth` function is a multi-threaded implementation of the Aberth-Ehrlich method for finding
 * the roots of a polynomial.
 *
 * Aberth's method is a method for finding the roots of a polynomial that is
 * robust but requires complex arithmetic even if the polynomial is real. This
 * is because it starts with complex initial approximations.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of a
 * polynomial. Each element of the vector corresponds to a term in the polynomial, starting from the
 * highest degree term and ending with the constant term. For example, if the polynomial is `3x^2 +
 * 2x +
 * @param[in] zs `zs` is a vector of complex numbers representing the initial guesses for the roots
 * of the polynomial. The function will update these values iteratively to converge to the actual
 * roots.
 * @param[in] options The `options` parameter is an object of type `Options` that contains the
 * maximum number of iterations (`max_iters`) and the tolerance (`tolerance`). These options control
 * the convergence criteria for the Aberth-Ehrlich method.
 *
 * @return The `aberth` function returns a `std::pair<unsigned int, bool>`. The first element of the
 * pair represents the number of iterations performed, and the second element represents whether the
 * method converged to a solution within the specified tolerance.
 */
extern auto aberth_autocorr_mt(const std::vector<double>& coeffs,
                               std::vector<std::complex<double>>& zs, const Options& options)
    -> std::pair<unsigned int, bool>;

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
