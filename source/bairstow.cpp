#include <cmath>                   // for abs, acos, cos, pow
#include <cstddef>                 // for size_t
#include <ginger/bairstow.hpp>     // for operator-, Vector2Ref
#include <ginger/rootfinding.hpp>  // for Vec2, delta, Options, horner_eval
#include <ginger/vector2.hpp>      // for operator-, Vector2
#include <iterator>                // for back_inserter
#include <utility>                 // for pair
#include <vector>                  // for vector, vector<>::reference, __v...

#ifndef M_PI
#    define M_PI 3.14159265358979323846264338327950288
#endif


// static double dummy_x = 0.0;
// static double dummy_y = 0.0;
// static const Vec2Ref dummy_ref(dummy_x, dummy_y); 

/**
 * The function `horner` implements the Horner's method for evaluating a
 * polynomial at a given point.
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
auto horner_ref(std::vector<double> &coeffs, std::vector<Vec2Ref> &vcoeffs, size_t degree,
                 const Vec2 &vr) -> Vec2Ref {
    auto itr0 = coeffs.begin();
    auto itr1 = std::next(vcoeffs.begin());
    for (auto i = 0U; i != degree - 1; ++i, ++itr0, ++itr1) {
        *itr1 += vr * (*itr0);
    }
    return vcoeffs[degree - 1];
}

/**
 * @brief Bairstow's method (even degree only)
 *
 * The `bairstow_even` function implements Bairstow's method for finding the roots of a real
 * polynomial with an even degree using multi-threading.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of the
 * polynomial. Each element of the vector corresponds to the coefficient of a term in the
 * polynomial, starting from the highest degree term and ending with the constant term. For example,
 * if the polynomial is `3x^2 +
 * 2
 * @param[in, out] vrs `vrs` is a vector of iterates, which represents the initial guesses for the
 * roots of the polynomial. The Bairstow's method will update these iterates iteratively until the
 * desired tolerance is reached or the maximum number of iterations is reached.
 * @param[in] options The `options` parameter is an object of type `Options` which contains the
 * maximum number of iterations (`max_iters`) and the tolerance (`tolerance`). These options are
 * used to control the convergence criteria for the Bairstow's method.
 *
 * @return The function `pbairstow_even` returns a `std::pair<unsigned int, bool>`. The first
 * element of the pair represents the number of iterations performed, and the second element
 * represents whether the method converged to a solution within the specified tolerance.
 */
auto bairstow(const std::vector<double> &coeffs, Vec2 &vr,
                   const Options &options = Options()) -> std::pair<unsigned int, bool> {
    auto coeffs1 = coeffs;
    const auto degree = coeffs1.size() - 1;  // degree, assume even
    std::vector<Vec2Ref> vcoeffs1;
    for (auto i = 0U; i < degree; ++i) {
        vcoeffs1.emplace_back(Vec2Ref{coeffs1[i], coeffs1[i+1]});
    }

    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        std::copy(coeffs.begin(), coeffs.end(), std::back_inserter(coeffs1));
        auto vA = horner_ref(coeffs1, vcoeffs1, degree, vr);
        auto vA1 = horner_ref(coeffs1, vcoeffs1, degree - 2, vr);
        const auto tolerance = std::max(std::abs(vA.x()), std::abs(vA.y()));
        vr -= delta_ref(vA, vr, vA1);
        if (tolerance < options.tolerance) {
            return {niter, true};
        }
    }
    return {options.max_iters, false};
}
