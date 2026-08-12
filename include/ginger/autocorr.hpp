/**
 * @file autocorr.hpp
 * @brief Auto-correlation polynomial root-finding (palindromic polynomials)
 */

#pragma once

#include <utility>
#include <vector>

#include "rootfinding.hpp"

/**
 * @brief Initial guess for the parallel Bairstow method (specific for auto-correlation function)
 *
 * The function `initial_autocorr` calculates the initial quadratic factors for finding the roots of
 * a palindromic (auto-correlation) polynomial.
 *
 * @f[
 *     x^2 - r_k x - q_k, \quad r_k = 2\cos(2\pi\phi_2(k)), \quad q_k = -1, \quad k =
 * 0,\dots,\lfloor n/2\rfloor
 * @f]
 *
 * @dot
 *   digraph autocorr_init {
 *     rankdir=LR; bgcolor="transparent";
 *     node [shape=box, style=filled, fillcolor="#d4e6f1"];
 *     coeffs [label="Polynomial\ncoeffs [a0..an]", fillcolor="#a9cce3"];
 *     vdc [label="VdC(2)\ncircle points"];
 *     init [label="Initial factors\nr_k, q_k\non unit circle", fillcolor="#f9e79f"];
 *     autocorr [label="Palindromic\nsymmetry", fillcolor="#d5f5e3"];
 *     coeffs -> init;
 *     vdc -> init;
 *     init -> autocorr;
 *   }
 * @enddot
 *
 * @param[in] coeffs The parameter `coeffs` is a vector of doubles.
 *
 * @return The function `initial_autocorr` returns a vector of `Vec2` objects.
 */
extern auto initial_autocorr(const std::vector<double>& coeffs) -> std::vector<Vec2>;

/**
 * @brief Single-threading Bairstow's method (specific for auto-correlation function)
 *
 * The function `pbairstow_autocorr_st` implements Bairstow's method for finding the roots of a
 * palindromic (auto-correlation) polynomial using a single thread.
 *
 * Each iterate finds a quadratic factor @f$ x^2 - r_i x - q_i @f$ that respects the palindromic
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
extern auto pbairstow_autocorr_st(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                                  const Options& options) -> std::pair<unsigned int, bool>;

inline auto pbairstow_autocorr(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                               const Options& options) -> std::pair<unsigned int, bool> {
    return pbairstow_autocorr_st(coeffs, vrs, options);
}

/**
 * @brief Extract autocorrelation quadratic factor
 *
 * Converts a quadratic factor of the form @f$ x^2 - r x - q @f$
 * to the autocorrelation form @f$ (-1/q) + (r/q) x + x^2 @f$
 * where the roots appear in reciprocal pairs.
 *
 * Quadratic factor relations:
 * @f[
 *     x^2 - r x - q = (x - a_1)(x - a_2) = x^2 - (a_1 + a_2)x + a_1 a_2
 * @f]
 * @f[
 *     x^2 + r x + t = (x + a_1)(x + a_2) = x^2 + (a_1 + a_2)x + a_1 a_2
 * @f]
 *
 * For autocorrelation polynomials, the transformed factor has roots @f$ a_i @f$ and @f$ 1/a_i @f$.
 *
 * @param[in,out] vr Quadratic factor Vec2 (r, q) modified in-place
 */
extern void extract_autocorr(Vec2& vr);

/**
 * @brief Reconstruct a monic polynomial from its autocorrelation quadratic factors
 *
 * Auto-correlation (palindromic) polynomials have roots in reciprocal pairs.
 * Each quadratic factor @f$ x^2 - r_i x - q_i @f$ found by pbairstow_autocorr carries 2 roots.
 * This function adds the reciprocal of each root, then reconstructs the full
 * monic polynomial with Leja ordering for numerical accuracy.
 *
 * @f[
 *     P(x) = \prod_{i=1}^{n/2} (x - a_i)(x - a_i^{-1}) = x^{2n} + a_{2n-1}x^{2n-1} + \cdots + a_0
 * @f]
 *
 * @param[in] vrs Quadratic factors from pbairstow_autocorr
 * @return std::vector<double> Monic polynomial coefficients (highest degree first)
 */
extern auto poly_from_autocorr_factors(const std::vector<Vec2>& vrs) -> std::vector<double>;
