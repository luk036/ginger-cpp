/**
 * @file bairstow.hpp
 * @brief Bairstow's root-finding method with reference-based vectors
 */

#pragma once

#include <utility>
#include <vector>

#include "matrix2.hpp"
#include "vector2.hpp"
// #include "vector2_ref.hpp"

using Vec2 = ginger::Vector2<double>;
using Mat2 = ginger::Matrix2<Vec2>;
using Vec2Ref = ginger::Vector2<double&>;

class Options;

/**
 * @brief Horner's rule
 *
 * Horner's rule is a method for evaluating a polynomial of degree degree at a given
 * point x. It involves rewriting the polynomial as a nested multiplication and
 * addition of the form:
 *
 *  P(x) = a_0 + x(a_1 + x(a_2 + ... + x(a_{degree-1} + x(a_n))...))
 *
 * This form allows for efficient evaluation of the polynomial at a given point
 * x using only degree multiplications and degree additions. Horner's rule is commonly
 * used in numerical methods for polynomial evaluation and interpolation.
 *
 * @param[in, out] coeffs coeffs is a reference to a vector of doubles. It is used to
 * store the coefficients of a polynomial.
 * @param[in, out] vcoeffs vcoeffs is a reference to a vector of Vec2Ref objects.
 * @param[in] degree The parameter `degree` represents the size of the vector `coeffs`. It
 * indicates the number of elements in the vector `coeffs`.
 * @param[in] vr vr is a Vec2 object, which represents a 2D vector. It has two
 * components, vr.x() and vr.y(), which are used in the calculations inside the
 * horner function.
 *
 * @return a Vec2Ref object.
 */
extern auto horner_ref(std::vector<double>& coeffs, std::vector<Vec2Ref>& vcoeffs,
                       std::size_t degree, const Vec2& vr) -> Vec2Ref;

/**
 * @brief Bairstow's method
 *
 * The `bairstow` function implements Bairstow's method for finding the roots of a real
 * polynomial.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of the
 * polynomial. Each element of the vector corresponds to the coefficient of a term in the
 * polynomial, starting from the highest degree term and ending with the constant term. For example,
 * if the polynomial is `3x^2 + 2x + 1`, the coefficients vector would be `{3, 2, 1}`.
 * @param[in,out] vr `vr` is a vector of iterates, which represents the initial guesses for the
 * roots of the polynomial. Bairstow's method will update these iterates iteratively until the
 * desired tolerance is reached or the maximum number of iterations is reached.
 * @param[in] options The `options` parameter is an object of type `Options` which contains the
 * maximum number of iterations (`max_iters`) and the tolerance (`tolerance`). These options are
 * used to control the convergence criteria for Bairstow's method.
 *
 * @return The function `bairstow` returns a `std::pair<unsigned int, bool>`. The first
 * element of the pair represents the number of iterations performed, and the second element
 * represents whether the method converged to a solution within the specified tolerance.
 */
extern auto bairstow(const std::vector<double>& coeffs, Vec2& vr, const Options& options)
    -> std::pair<unsigned int, bool>;

/**
 * @brief Create adjoint matrix from two vectors (reference-based)
 *
 * @param[in] vr Vector r (constant reference)
 * @param[in] vp Vector p (reference wrapper)
 * @return Mat2 The adjoint 2x2 matrix
 */
inline auto makeadjoint_ref(const Vec2& vr, const Vec2Ref& vp) -> Mat2 {
    auto p = vp.x();
    auto s = vp.y();
    return {Vec2{s, -p}, Vec2{-p * vr.y(), p * vr.x() + s}};
}

/**
 * @brief Calculate Newton correction delta (reference-based)
 *
 * Uses the adjoint matrix to compute the correction step in Bairstow's method.
 *
 * @param[in] vA Current remainder vector (reference wrapper)
 * @param[in] vr Direction vector r
 * @param[in] vp Vector p (reference wrapper)
 * @return Vec2 The correction delta
 */
inline auto delta_ref(const Vec2Ref& vA, const Vec2& vr, const Vec2Ref& vp) -> Vec2 {
    const auto mp = makeadjoint_ref(vr, vp);  // 2 mul's
    return mp.mdot(vA) / mp.det();            // 6 mul's + 2 div's
}
