#include <ginger/ThreadPool.h>  // for ThreadPool

#include <cmath>    // for abs, acos, cos, pow
#include <cstddef>  // for size_t
#include <ginger/config.hpp>
#include <ginger/robin.hpp>        // for Robin
#include <ginger/rootfinding.hpp>  // for Vec2, delta, Options, horner_eval
#include <ginger/vector2.hpp>      // for operator-, Vector2
#include <utility>  // for pair
#include <vector>   // for vector, vector<>::reference, __v...

#ifndef M_PI
#    define M_PI 3.14159265358979323846264338327950288
#endif

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
 *
 * ```svgbob
 * Horner's method for polynomial evaluation:
 *
 * For polynomial P(x) = a_n*x^n + a_(n-1)*x^(n-1) + ... + a_1*x + a_0
 *
 * The algorithm proceeds as:
 * t = 0
 * for i = 0 to n-1:
 *   t = t + coeffs[i] * vr.x()  // accumulate x coefficient
 *   t = t + coeffs[i] * vr.y()  // accumulate y coefficient (for quadratic factor)
 *
 * coeffs[0]     coeffs[1]     coeffs[2]               coeffs[n-2]    coeffs[n-1]
 * +-----------> +-----------> +-----------> + ... + -----------> + --------->
 * |             |             |                      |            |
 * |    vr.x()   v    vr.x()   v                      v    vr.x()   v    vr.x()
 * +---> [x] --> +---> [x] --> +---> [...] --> ... +---> [x] --> +---> [x] -->
 * |             |             |                      |            |
 * |    vr.y()   v    vr.y()   v                      v    vr.y()   v    vr.y()
 * +---> [x] --> +---> [x] --> +---> [...] --> ... +---> [x] --> +---> [x] -->
 * |             |             |                      |            |
 * +-------------+-------------+----------------------+------------+-----------> Vec2 result
 *
 * Where each step modifies the coefficients in place
 * ```
 */
auto horner(std::vector<double> &coeffs1, size_t degree, const Vec2 &vr) -> Vec2 {
    auto itr0 = coeffs1.begin();
    auto itr1 = std::next(itr0);
    auto itr2 = std::next(itr1);
    for (auto i = 0U; i != degree - 1; ++i, ++itr0, ++itr1, ++itr2) {
        *itr1 += *itr0 * vr.x();
        *itr2 += *itr0 * vr.y();
        // coeffs1[i + 1] += coeffs1[i] * vr.x();
        // coeffs1[i + 2] += coeffs1[i] * vr.y();
    }
    return Vec2{coeffs1[degree - 1], coeffs1[degree]};
}

/**
 * The function `suppress` calculates and updates the values of `vA` and `vA1`
 * based on the given input vectors `vri` and `vrj`.
 *
 * @param[in, out] vA A reference to a Vec2 object representing vector A.
 * @param[in, out] vA1 vA1 is a reference to a Vec2 object.
 * @param[in] vri A vector representing the position of point i.
 * @param[in] vrj The parameter `vrj` represents a `Vec2` object.
 *
 * ```svgbob
 * Suppression step in root-finding:
 *
 * Given two root approximations vri and vrj, suppress the influence of vrj on vri
 *
 *    vri = (r_i, q_i)  <-- current root approximation
 *         |
 *         |  vp = vri - vrj  (difference vector)
 *         v
 *      +--------+    Calculate adjoint matrix and determinant
 *      | adjoint|----> m_adjoint = [[s, -p], [-p*vri.y, p*vri.x+s]]
 *      | matrix |                 where p=vp.x, s=vp.y
 *      +--------+ 
 *         |
 *         v
 *      +--------+    Apply matrix transformation to update vA, vA1
 *      | update |----> vA, vA1 = transformed values
 *      | vA, vA1|
 *      +--------+
 *         |
 *         v
 *    Updated vA, vA1 without influence from vrj
 *
 * This process removes the contribution of root vrj from the evaluation at vri
 * ```
 */
auto suppress(Vec2 &vA, Vec2 &vA1, const Vec2 &vri, const Vec2 &vrj) -> void {
    const auto vp = vri - vrj;
    const auto p = vp.x();
    const auto s = vp.y();
    const auto m_adjoint = Mat2{Vec2{s, -p}, Vec2{-p * vri.y(), p * vri.x() + s}};
    const auto e = m_adjoint.det();
    const auto va = m_adjoint.mdot(vA);
    const auto vd = vA1 * e - va;
    const auto vc = Vec2{vd.x(), vd.y() - va.x() * p};
    vA = va * e;
    vA1 = m_adjoint.mdot(vc);
}

/**
 * The function `suppress2` calculates and updates the values of `vA` and `vA1`
 * based on the given input vectors `vri` and `vrj`.
 *
 * @param[in, out] vA A reference to a Vec2 object representing vector A.
 * @param[in, out] vA1 vA1 is a reference to a Vec2 object.
 * @param[in] vri A vector representing the position of point i.
 * @param[in] vrj The parameter `vrj` represents a `Vec2` object.
 *
 * ```svgbob
 * Alternative suppression step in root-finding:
 *
 * This version uses vrj in the adjoint calculation instead of vri
 *
 *    vri = (r_i, q_i), vrj = (r_j, q_j)  <-- root approximations
 *         |
 *         |  vp = vri - vrj  (difference vector)
 *         v
 *      +--------+    Calculate adjoint matrix using vrj instead of vri
 *      | adjoint|----> m_adjoint = [[s, -p], [-p*vrj.y, p*vrj.x+s]]
 *      | matrix |                 where p=vp.x, s=vp.y
 *      +--------+ 
 *         |
 *         v
 *      +--------+    Scale and adjust vA, vA1 values
 *      | update |----> vA *= e, vA1 *= e, vA1 -= m_adjoint.mdot(vA)
 *      | vA, vA1|
 *      +--------+
 *         |
 *         v
 *    Updated vA, vA1 with alternative suppression method
 *
 * Specialized variant for specific root-finding scenarios
 * ```
 */
auto suppress2(Vec2 &vA, Vec2 &vA1, const Vec2 &vri, const Vec2 &vrj) -> void {
    const auto vp = vri - vrj;
    auto p = vp.x();
    auto s = vp.y();
    const auto m_adjoint = Mat2{Vec2{s, -p}, Vec2{-p * vrj.y(), p * vrj.x() + s}};
    const auto e = m_adjoint.det();
    vA1 *= e;  // e may tend to zero
    vA1 -= m_adjoint.mdot(vA);
    vA *= e;
}

/**
 * @brief Initial guess for the parallel Bairstow method
 *
 * The `initial_guess` function calculates the initial values for the parallel Bairstow method for
 * finding the roots of a real polynomial.
 *
 * @param[in] coeffs coeffs is a vector of doubles that represents the coefficients of a polynomial.
 *
 * @return The function `initial_guess` returns a vector of `Vec2` objects.
 *
 * ```svgbob
 * Initial guess calculation for Bairstow's method:
 *
 *                    center
 *                      *
 *                     /|\
 *                    / | \ radius
 *                   /  |  \
 *    (-2*radius,-m) *   |   * (2*radius,-m)
 *                   \  |  /
 *                    \ | /
 *                     \|/
 *                      *
 *                     (0,-m)
 *
 * Where:
 * - center = -coeffs[1] / (degree * coeffs[0])
 * - radius = |P(center)|^(1/degree)
 * - m = center^2 + radius^2
 * - Points placed at (2*(center + radius*cos(θ)), -(m + 2*center*radius*cos(θ)))
 * - θ values: π*i/degree for odd i from 1 to degree-1
 * ```
 */
auto initial_guess(std::vector<double> coeffs) -> std::vector<Vec2> {
    auto degree = coeffs.size() - 1;
    const auto center = -coeffs[1] / (double(degree) * coeffs[0]);
    const auto poly_c = horner_eval(std::move(coeffs), degree, center);
    const auto radius = std::pow(std::abs(poly_c), 1.0 / double(degree));
    degree /= 2;
    degree *= 2;  // make even
    const auto k = M_PI / double(degree);
    const auto m = center * center + radius * radius;
    auto vr0s = std::vector<Vec2>{};
    for (auto i = 1U; i < degree; i += 2) {
        const auto temp = radius * std::cos(k * i);
        auto r0 = 2 * (center + temp);
        auto t0 = -(m + 2 * center * temp);
        vr0s.emplace_back(Vec2{r0, t0});
    }
    return vr0s;
}

/**
 * @brief Multi-threading Bairstow's method (even degree only)
 *
 * The `pbairstow_even` function implements Bairstow's method for finding the roots of a real
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
 *
 * ```svgbob
 * Parallel Bairstow's method iterative process:
 *
 * For each iterate vr_i, the process is:
 *
 *  coeffs +--------+ P(vr_i) 
 *         | horner |------>
 *         |        | P'(vr_i) 
 *         +--------+------>
 *              |
 *              v
 *         +--------+     For each other iterate vr_j (j ≠ i):
 *         |suppress|----> Remove vr_j's influence on vr_i
 *         |        |     (apply suppress transformation)
 *         +--------+
 *              |
 *              v
 *         +--------+ 
 *         | update | vr_i^(k+1) = vr_i^(k) - delta(P(vr_i), vr_i, P'(vr_i))
 *         | vr_i   | considering all other roots vr_j (j ≠ i)
 *         +--------+
 *              |
 *              v
 *         vr_i^(k+1)
 *
 * All iterates updated in parallel using thread pool
 * Convergence check across all iterates simultaneously
 * ```
 */
auto pbairstow_even(const std::vector<double> &coeffs, std::vector<Vec2> &vrs,
                    const Options &options = Options()) -> std::pair<unsigned int, bool> {
    ThreadPool pool(std::thread::hardware_concurrency());

    const auto M = vrs.size();
    const auto rr = fun::Robin<size_t>(M);

    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        std::vector<std::future<double>> results;

        for (auto i = 0U; i != M; ++i) {
            results.emplace_back(pool.enqueue([&coeffs, &vrs, &rr, i]() {
                const auto degree = coeffs.size() - 1;  // degree, assume even
                const auto &vri = vrs[i];
                auto coeffs1 = coeffs;
                auto vA = horner(coeffs1, degree, vri);
                auto vA1 = horner(coeffs1, degree - 2, vri);
                const auto tol_i = std::max(std::abs(vA.x()), std::abs(vA.y()));
                for (auto j : rr.exclude(i)) {
                    const auto vrj = vrs[j];  // make a copy, don't reference!
                    suppress(vA, vA1, vri, vrj);
                }
                vrs[i] -= delta(vA, vri, vA1);  // Gauss-Seidel fashion
                return tol_i;
            }));
        }
        for (auto &&result : results) {
            auto &&res = result.get();
            if (tolerance < res) {
                tolerance = res;
            }
        }
        if (tolerance < options.tolerance) {
            return {niter, true};
        }
    }
    return {options.max_iters, false};
}
