#pragma once

// import numpy as np
#include <utility>
#include <vector>

#include "matrix2.hpp"
#include "vector2.hpp"

using Vec2 = ginger::Vector2<double>;
using Mat2 = ginger::Matrix2<Vec2>;

class Options;

/**
 * The function calculates the initial autocorrelation values (specific for
 * auto-correlation function)
 *
 * @param[in] coeffs The parameter `coeffs` is a vector of doubles.
 *
 * @return The function `initial_autocorr` returns a vector of `Vec2` objects.
 */
extern auto initial_autocorr(const std::vector<double>& coeffs) -> std::vector<Vec2>;

/**
 * @brief Multi-threading Bairstow's method (specific for auto-correlation
 * function)
 *
 * The function `pbairstow_autocorr` is implementing the Bairstow's method for
 * finding the roots of a real polynomial (specific for auto-correlation
 * function)
 *
 * @param[in] coeffs polynomial
 * @param[in,out] vrs vector of iterates
 * @param[in] options maximum iterations and tolorance
 * @return std::pair<unsigned int, bool>
 */
extern auto pbairstow_autocorr(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                               const Options& options) -> std::pair<unsigned int, bool>;

/**
 * The function extracts the autocorrelation values from a given vector.
 *
 *   x^2 - r*x - q  or (-1/q) + (r/q) * x + x^2
 *   (x - a1)(x - a2) = x^2 - (a1 + a2) x + a1 * a2
 *
 *   x^2 + r*x + t or x^2 + (r/t) * x + (1/t)
 *   (x + a1)(x + a2) = x^2 + (a1 + a2) x + a1 * a2
 *
 * @param[in,out] vr The parameter `vr` is of type `Vec2`, which is a custom
 * class representing a 2D vector. It contains two components, `x` and `y`,
 * which are accessed using the `x()` and `y()` member functions respectively.
 */
extern void extract_autocorr(Vec2& vr);

/**
 * @brief Reconstruct a monic polynomial from its autocorrelation quadratic factors
 *
 * Auto-correlation (palindromic) polynomials have roots in reciprocal pairs.
 * Each quadratic factor x^2 - r*x - q found by pbairstow_autocorr carries 2 roots.
 * This function adds the reciprocal of each root, then reconstructs the full
 * monic polynomial with Leja ordering for numerical accuracy.
 *
 * @param[in] vrs Quadratic factors from pbairstow_autocorr
 * @return std::vector<double> Monic polynomial coefficients (highest degree first)
 */
extern auto poly_from_autocorr_factors(const std::vector<Vec2>& vrs) -> std::vector<double>;
