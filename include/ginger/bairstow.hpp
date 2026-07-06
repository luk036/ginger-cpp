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
 * @brief Horner's rule (reference-based Vec2 version)
 *
 * Horner's rule evaluates a polynomial at a given point \f$x\f$.
 * It rewrites the polynomial as nested multiplication:
 * @f[
 *     P(x) = a_0 + x(a_1 + x(a_2 + \cdots + x(a_{n-1} + x a_n)\cdots))
 * @f]
 * This allows evaluation using \f$n\f$ multiplications and \f$n\f$ additions.
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
 * Bairstow's method finds quadratic factors of the form \f$x^2 - rx - q\f$ of a real polynomial
 * using Newton's method in 2D. At each iteration, the correction is:
 * @f[
 *     \begin{bmatrix} \Delta r \\ \Delta q \end{bmatrix} = -J^{-1} \begin{bmatrix} P(r,q) \\ Q(r,q)
 * \end{bmatrix}
 * @f]
 * where \f$P\f$ and \f$Q\f$ are the remainders of synthetic division by \f$x^2 - rx - q\f$.
 *
 * @dot
 *   digraph bairstow_iter {
 *     rankdir=LR; bgcolor="transparent";
 *     node [shape=box, style=filled, fillcolor="#d4e6f1"];
 *     init [label="Initial\n(r, q)", fillcolor="#a9cce3"];
 *     horner [label="Synthetic\ndivision"];
 *     jacobian [label="Jacobian\n+ adjoint"];
 *     newton [label="Newton step\nDelta r, Delta q"];
 *     check [label="|P|,|Q| < tol?", shape=diamond, fillcolor="#f9e79f"];
 *     done [label="Quadratic\nfactor found!", fillcolor="#7fb3d8"];
 *     init -> horner -> jacobian -> newton -> check;
 *     check -> horner [label="No", style=dashed, color="#e74c3c"];
 *     check -> done [label="Yes", color="#27ae60"];
 *   }
 * @enddot
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
 * Computes the adjugate matrix of the Jacobian in Bairstow's method:
 * @f[
 *     \operatorname{adj}(J) = \begin{bmatrix} s & -p \cdot r_y \\ -p & p \cdot r_x + s
 * \end{bmatrix}
 * @f]
 * where @f$ (p, s) = vp @f$ and @f$ (r_x, r_y) = vr @f$.
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
 * Uses the adjoint matrix to compute the correction step in Bairstow's method:
 * @f[
 *     \begin{bmatrix} \Delta r \\ \Delta q \end{bmatrix} = -\frac{\operatorname{adj}(J)}{\det(J)}
 * \, vA
 * @f]
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
