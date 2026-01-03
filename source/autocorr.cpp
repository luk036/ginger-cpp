#include <ginger/ThreadPool.h>  // for ThreadPool

#include <cmath>    // for abs, acos, cos, pow
#include <cstddef>  // for size_t
#include <future>   // for future
#include <ginger/config.hpp>
#include <ginger/robin.hpp>        // for Robin
#include <ginger/rootfinding.hpp>  // for Vec2, delta, Options, horner_eval
#include <ginger/vector2.hpp>      // for operator-, Vector2
#include <thread>                  // for thread
#include <utility>                 // for pair
#include <vector>                  // for vector, vector<>::reference, __v...

#ifndef M_PI
#    define M_PI 3.14159265358979323846264338327950288
#endif

/**
 * The function calculates the initial autocorrelation values (specific for
 * auto-correlation function)
 *
 * @param[in] coeffs The parameter `coeffs` is a vector of doubles.
 *
 * @return The function `initial_autocorr` returns a vector of `Vec2` objects.
 *
 * ```svgbob
 * For auto-correlation functions:
 *
 * Initial values distributed as:
 *
 *        real
 *         |
 *         |
 *   *-----+-----*----->
 *   |     |     |
 *   |     |     | radius^2
 *   |     |     |
 *   *-----+-----*
 *         |
 *         v
 *        imag
 *
 * Where radius = |coeffs[degree]|^(1/degree)
 * Points placed at specific positions for auto-correlation property
 * ```
 */
auto initial_autocorr(const std::vector<double> &coeffs) -> std::vector<Vec2> {
    auto degree = coeffs.size() - 1;
    const auto radius = std::pow(std::abs(coeffs[degree]), 1.0 / double(degree));

    degree /= 2;
    const auto k = M_PI / double(degree);
    const auto m = radius * radius;
    auto vr0s = std::vector<Vec2>{};
    for (auto i = 1U; i < degree; i += 2) {
        vr0s.emplace_back(Vec2{2 * radius * std::cos(k * double(i)), -m});
    }
    return vr0s;
}

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
 *
 * ```svgbob
 * Bairstow's method iterative process for auto-correlation:
 *
 * For each iterate vr_i, the process is:
 *
 *  coeffs +--------+ P(vr_i)
 *         | horner |------>
 *         |        | P'(vr_i)
 *         +--------+------>
 *              |
 *              v
 *         +--------+
 *         | update | vr_i^(k+1) = vr_i^(k) - [P(vr_i^(k)) / P'(vr_i^(k))]
 *         | vr_i   | considering other roots vr_j (j ≠ i) and special auto-correlation terms
 *         +--------+
 *              |
 *              v
 *         vr_i^(k+1)
 *
 * Parallel computation across all iterates vr_0, vr_1, ..., vr_n
 * Special handling for auto-correlation property: process both vr_j and (1/vr_j)
 * ```
 */
auto pbairstow_autocorr(const std::vector<double> &coeffs, std::vector<Vec2> &vrs,
                        const Options &options = Options()) -> std::pair<unsigned int, bool> {
    ThreadPool pool(std::thread::hardware_concurrency());

    const auto num_roots = vrs.size();
    const auto rr = fun::Robin<size_t>(num_roots);

    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        std::vector<std::future<double>> results;
        for (auto idx = 0U; idx != num_roots; ++idx) {
            results.emplace_back(pool.enqueue([&, idx]() {
                auto coeffs1 = coeffs;
                const auto degree = coeffs.size() - 1;  // degree, assume even
                const auto &vri = vrs[idx];
                auto vA = horner(coeffs1, degree, vri);
                const auto tol_i = std::max(std::abs(vA.x()), std::abs(vA.y()));
                auto vA1 = horner(coeffs1, degree - 2, vri);
                for (auto jdx : rr.exclude(idx)) {  // exclude idx
                    const auto vrj = vrs[jdx];    // make a copy, don't reference!
                    suppress(vA, vA1, vri, vrj);
                    const auto vrjn = ginger::Vector2<double>(-vrj.x(), 1.0) / vrj.y();
                    suppress(vA, vA1, vri, vrjn);
                }
                const auto vrin = ginger::Vector2<double>(-vri.x(), 1.0) / vri.y();
                suppress(vA, vA1, vri, vrin);

                vrs[idx] -= delta(vA, vri, std::move(vA1));  // Gauss-Seidel fashion
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
 *
 * ```svgbob
 * Root extraction process:
 *
 * Input: vr = (r, q) representing x^2 - r*x - q
 *
 * Calculate:
 *   h = r/2
 *   d = h^2 + q
 *
 *       d >= 0?  ----->  yes  --------->   calculate real roots
 *          |                                  a1 = h + sign(h)*sqrt(d)
 *         no                                 a2 = -q / a1
 *          |                              check if |a1| > 1 or |a2| > 1
 *          v                              if so, replace with 1/a1 or 1/a2
 *    complex roots                        update vr = (a1+a2, -a1*a2)
 *    check if |q| > 1
 *    if so, transform
 *    vr = (-r/q, 1/q)
 *
 * Output: updated vr
 * ```
 */
void extract_autocorr(Vec2 &vr) {
    const auto &r = vr.x();
    const auto &q = vr.y();
    const auto hr = r / 2.0;
    const auto d = hr * hr + q;
    if (d < 0.0) {  // complex conjugate root
        if (q < -1.0) {
            vr = Vec2{-r, 1.0} / q;
        }
        // else no need to change
    } else {  // two real roots
        auto a1 = hr + (hr >= 0.0 ? sqrt(d) : -sqrt(d));
        auto a2 = -q / a1;
        if (std::abs(a1) > 1.0) {
            if (std::abs(a2) > 1.0) {
                a2 = 1.0 / a2;
            }
            a1 = 1.0 / a1;
            vr = Vec2{a1 + a2, -a1 * a2};
        } else if (std::abs(a2) > 1.0) {
            a2 = 1.0 / a2;
            vr = Vec2{a1 + a2, -a1 * a2};
        }
        // else no need to change
    }
}
