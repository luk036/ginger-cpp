/**
 * @file rootfinding.hpp
 * @brief Parallel Bairstow root-finding methods for real polynomials
 */

#pragma once

#include <utility>
#include <vector>

#include "matrix2.hpp"
#include "vector2.hpp"

using Vec2 = ginger::Vector2<double>;
using Mat2 = ginger::Matrix2<Vec2>;

/**
 * @brief Options
 *
 * The code snippet defines a class called `Options` that represents the options for a specific
 * algorithm or function. It has two public member variables: `max_iters` and `tolerance`.
 */
class Options;

/**
 * @brief Initial guess for the parallel Bairstow method
 *
 * The `initial_guess` function calculates the initial values for the parallel Bairstow method for
 * finding the roots of a real polynomial.
 *
 * @param[in] coeffs coeffs is a vector of doubles that represents the coefficients of a polynomial.
 * The vector is passed by value.
 *
 * @return The function `initial_guess` returns a vector of `Vec2` objects.
 */
extern auto initial_guess(std::vector<double> coeffs) -> std::vector<Vec2>;

/**
 * @brief Multi-threading Bairstow's method (even degree only)
 *
 * The `pbairstow_even` function implements Bairstow's method for finding the roots of a real
 * polynomial with an even degree using multi-threading.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of the
 * polynomial. Each element of the vector corresponds to the coefficient of a term in the
 * polynomial, starting from the highest degree term and ending with the constant term. For example,
 * if the polynomial is `3x^2 + 2x + 1`, the coefficients vector would be `{3, 2, 1}`.
 * @param[in,out] vrs `vrs` is a vector of iterates, which represents the initial guesses for the
 * roots of the polynomial. Bairstow's method will update these iterates iteratively until the
 * desired tolerance is reached or the maximum number of iterations is reached.
 * @param[in] options The `options` parameter is an object of type `Options` which contains the
 * maximum number of iterations (`max_iters`) and the tolerance (`tolerance`). These options are
 * used to control the convergence criteria for Bairstow's method.
 *
 * @return The function `pbairstow_even` returns a `std::pair<unsigned int, bool>`. The first
 * element of the pair represents the number of iterations performed, and the second element
 * represents whether the method converged to a solution within the specified tolerance.
 */
extern auto pbairstow_even(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                           const Options& options) -> std::pair<unsigned int, bool>;

/**
 * @brief Horner's rule
 *
 * Horner's rule is a method for evaluating a polynomial at a given point \f$x\f$.
 * It rewrites the polynomial as nested multiplication:
 * @f[
 *     P(x) = a_0 + x(a_1 + x(a_2 + \cdots + x(a_{n-1} + x a_n)\cdots))
 * @f]
 * This allows evaluation using \f$n\f$ multiplications and \f$n\f$ additions.
 *
 * @param[in, out] coeffs1 coeffs1 is a reference to a vector of doubles. It is used to
 * store the coefficients of a polynomial.
 * @param[in] degree The parameter `degree` represents the size of the vector `coeffs1`. It
 * indicates the number of elements in the vector `coeffs1`.
 * @param[in] vr vr is a Vec2 object, which represents a 2D vector. It has two
 * components, vr.x() and vr.y(), which are used in the calculations inside the
 * horner function.
 *
 * @return a Vec2 object.
 */
extern auto horner(std::vector<double>& coeffs1, std::size_t degree, const Vec2& vr) -> Vec2;

/**
 * @brief Zero suppression step in Bairstow's method (variant 1)
 *
 * Finds coefficients of the linear remainder of a deflated polynomial
 * without explicitly constructing the deflated polynomial, avoiding
 * complex arithmetic within iterations.
 *
 * @param[in,out] vA First remainder coefficient
 * @param[in,out] vA1 Second remainder coefficient
 * @param[in] vri First known factor
 * @param[in] vrj Second known factor
 */
extern auto suppress(Vec2& vA, Vec2& vA1, const Vec2& vri, const Vec2& vrj) -> void;

/**
 * @brief Zero suppression step in Bairstow's method (variant 2)
 *
 * Alternative formulation of the zero suppression technique.
 *
 * @param[in,out] vA First remainder coefficient
 * @param[in,out] vA1 Second remainder coefficient
 * @param[in] vri First known factor
 * @param[in] vrj Second known factor
 */
extern auto suppress2(Vec2& vA, Vec2& vA1, const Vec2& vri, const Vec2& vrj) -> void;

/**
 * The function "makeadjoint" takes in a vector vr and a vector vp, and returns a 2x2 matrix where
 * the elements are calculated based on the values of vr and vp.
 *
 * @param[in] vr A constant reference to a Vec2 object, representing the vector vr.
 * @param[in] vp vp is a vector with two components, vp.x() and vp.y().
 *
 * @return a `Mat2` object.
 */
inline auto makeadjoint(const Vec2& vr, const Vec2& vp) -> Mat2 {
    auto&& p = vp.x();
    auto&& s = vp.y();
    return {Vec2{s, -p}, Vec2{-p * vr.y(), p * vr.x() + s}};
}

/**
 * The function calculates the delta value using the given parameters.
 *
 * @param[in] vA A vector of type Vec2.
 * @param[in] vr A vector representing the direction of rotation.
 * @param[in] vp The parameter `vp` is a `Vec2` object that is passed by rvalue reference.
 *
 * @return a Vec2 object.
 */
inline auto delta(const Vec2& vA, const Vec2& vr, const Vec2& vp) -> Vec2 {
    const auto mp = makeadjoint(vr, vp);  // 2 mul's
    return mp.mdot(vA) / mp.det();        // 6 mul's + 2 div's
}

/**
 * The function `horner_eval` evaluates a polynomial using Horner's method.
 *
 * Evaluates:
 * @f[
 *     P(z) = a_0 + a_1 z + a_2 z^2 + \cdots + a_n z^n
 * @f]
 * using the recurrence \f$b_0 = a_0, \; b_{k+1} = a_{k+1} + b_k z\f$, returning \f$b_n\f$.
 *
 * @param[in,out] coeffs1 A vector of coefficients for a polynomial, where the coefficient at index
 * i corresponds to the term with degree i.
 * @param[in] degree The degree parameter represents the degree of the polynomial. It indicates the
 * highest power of the variable in the polynomial equation.
 * @param[in] z The parameter `z` is a constant value that is used as the input to the polynomial
 * function being evaluated.
 *
 * @return a double value.
 */
inline auto horner_eval(std::vector<double> coeffs1, std::size_t degree, const double& z)
    -> double {
    for (auto idx = 0U; idx != degree; ++idx) {
        coeffs1[idx + 1] += coeffs1[idx] * z;
    }
    return coeffs1[degree];
}

/**
 * @brief Reconstruct a monic polynomial from its quadratic factors
 *
 * Given the quadratic factors found by Bairstow's method (each representing
 * \f$x^2 - r_i x - q_i\f$), multiply them together:
 * @f[
 *     P(x) = \prod_{i=1}^{n/2} (x^2 - r_i x - q_i) = x^n + a_{n-1}x^{n-1} + \cdots + a_0
 * @f]
 * To get the original polynomial, multiply the result by the original leading coefficient.
 *
 * @param[in] vrs Vector of quadratic factors, each as a Vec2 with x() = r, y() = q
 * @return std::vector<double> Monic polynomial coefficients (highest degree first)
 */
extern auto poly_from_quadratic_factors(const std::vector<Vec2>& vrs) -> std::vector<double>;
