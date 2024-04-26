#include <ginger/ThreadPool.h>  // for ThreadPool

#include <ginger/robin.hpp>        // for Robin
#include <ginger/rootfinding.hpp>  // for Vec2, delta, Options, horner_eval
#include <ginger/vector2.hpp>      // for operator-, Vector2
#include <cmath>                     // for abs, acos, cos, pow
#include <cstddef>                   // for size_t
// #include <functional>                // for __base
// #include <future>                    // for future
// #include <thread>                    // for thread
// #include <type_traits>               // for move
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
 * The function `suppress` calculates and updates the values of `vA` and `vA1`
 * based on the given input vectors `vri` and `vrj`.
 *
 * @param[in, out] vA A reference to a Vec2 object representing vector A.
 * @param[in, out] vA1 vA1 is a reference to a Vec2 object.
 * @param[in] vri A vector representing the position of point i.
 * @param[in] vrj The parameter `vrj` represents a `Vec2` object.
 */
auto suppress2(Vec2 &vA, Vec2 &vA1, const Vec2 &vri, const Vec2 &vrj) -> void {
    const auto vp = vri - vrj;
    auto p = vp.x();
    auto s = vp.y();
    const auto m_adjoint = Mat2{Vec2{s, -p}, Vec2{-p * vrj.y(), p * vrj.x() + s}};
    const auto e = m_adjoint.det();
    vA1 *= e; // e may tend to zero
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
 */
auto initial_guess(std::vector<double> coeffs) -> std::vector<Vec2> {
    auto N = coeffs.size() - 1;
    const auto c = -coeffs[1] / (double(N) * coeffs[0]);
    const auto Pc = horner_eval(std::move(coeffs), N, c);
    const auto re = std::pow(std::abs(Pc), 1.0 / double(N));
    N /= 2;
    N *= 2;  // make even
    const auto k = M_PI / double(N);
    const auto m = c * c + re * re;
    auto vr0s = std::vector<Vec2>{};
    for (auto i = 1U; i < N; i += 2) {
        const auto temp = re * std::cos(k * i);
        auto r0 = 2 * (c + temp);
        auto t0 = -(m + 2 * c * temp);
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
 * maximum number of iterations (`max_iters`) and the tolerance (`tolerance`). These options are used to
 * control the convergence criteria for the Bairstow's method.
 *
 * @return The function `pbairstow_even` returns a `std::pair<unsigned int, bool>`. The first
 * element of the pair represents the number of iterations performed, and the second element
 * represents whether the method converged to a solution within the specified tolerance.
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
