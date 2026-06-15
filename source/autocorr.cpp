#include <algorithm>
#include <cmath>              // for abs, acos, cos, pow
#include <complex>            // for complex
#include <cstddef>            // for size_t
#include <future>             // for future
#include <ginger/aberth.hpp>  // for poly_from_roots
#include <ginger/config.hpp>
#include <ginger/robin.hpp>        // for Robin
#include <ginger/rootfinding.hpp>  // for Vec2, delta, Options, horner_eval
#include <ginger/thread_pool.hpp>  // for thread_pool
#include <ginger/vector2.hpp>      // for operator-, Vector2
#include <thread>                  // for thread
#include <utility>                 // for pair
#include <vector>                  // for vector, vector<>::reference, __v...

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846264338327950288;
#endif

/**
 * The function calculates the initial autocorrelation values (specific for
 * auto-correlation function)
 *
 * @param[in] coeffs The parameter `coeffs` is a vector of doubles.
 *
 * @return The function `initial_autocorr` returns a vector of `Vec2` objects.
 *
 * @verbatim
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
 * @endverbatim
 */
auto initial_autocorr(const std::vector<double>& coeffs) -> std::vector<Vec2> {
    auto degree = coeffs.size() - 1;
    const auto radius = std::pow(std::abs(coeffs[degree]), 1.0 / static_cast<double>(degree));

    degree /= 2;
    const auto m = radius * radius;
    const auto num_points = degree / 2;
    auto vr0s = std::vector<Vec2>{};
    vr0s.reserve(num_points);
    for (auto i = 0U; i < num_points; ++i) {
        vr0s.emplace_back(2 * radius * cos_pi_vdc2(i), -m);
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
 * @verbatim
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
 * @endverbatim
 */
auto pbairstow_autocorr(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                        const Options& options = Options()) -> std::pair<unsigned int, bool> {
    auto& pool = ginger::get_thread_pool();

#if !defined(_MSC_VER) || !defined(_DEBUG)
    // MSVC Debug: thread_local at function scope has TLS guard init issues
    // when accessed from thread pool lambdas. Each task creates a local copy.
    thread_local std::vector<double> local_coeffs;
#endif

    const auto num_roots = vrs.size();
    const auto rr = fun::Robin<size_t>(num_roots);

    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        std::vector<std::future<double>> results;
        results.reserve(num_roots);
        for (auto idx = 0U; idx != num_roots; ++idx) {
            results.emplace_back(pool.enqueue([&coeffs, &vrs, &rr, idx]() {
                const auto degree = coeffs.size() - 1;
                const auto& vri = vrs[idx];
#if defined(_MSC_VER) && defined(_DEBUG)
                auto local_coeffs = coeffs;
#else
                local_coeffs = coeffs;
#endif
                auto vA = horner(local_coeffs, degree, vri);
                const auto tol_i = std::max(std::abs(vA.x()), std::abs(vA.y()));
                auto vA1 = horner(local_coeffs, degree - 2, vri);
                for (auto jdx : rr.exclude(idx)) {  // exclude idx
                    const auto vrj = vrs[jdx];      // make a copy, don't reference!
                    suppress(vA, vA1, vri, vrj);
                    const auto vrjn = ginger::Vector2<double>(-vrj.x(), 1.0) / vrj.y();
                    suppress(vA, vA1, vri, vrjn);
                }
                const auto vrin = ginger::Vector2<double>(-vri.x(), 1.0) / vri.y();
                suppress(vA, vA1, vri, vrin);

                vrs[idx] -= delta(vA, vri, vA1);  // Gauss-Seidel fashion
                return tol_i;
            }));
        }
        for (auto&& result : results) {
            auto&& res = result.get();
            tolerance = std::max(tolerance, res);
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
 * @verbatim
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
 * @endverbatim
 */
static auto roots_from_quadratic(const Vec2& vr)
    -> std::pair<std::complex<double>, std::complex<double>> {
    const auto r = vr.x();
    const auto q = vr.y();
    const auto disc = r * r + 4.0 * q;
    if (disc >= 0.0) {
        const auto sqrt_disc = std::sqrt(disc);
        return {{(r + sqrt_disc) / 2.0, 0.0}, {(r - sqrt_disc) / 2.0, 0.0}};
    }
    const auto sqrt_disc = std::sqrt(-disc);
    return {{r / 2.0, sqrt_disc / 2.0}, {r / 2.0, -sqrt_disc / 2.0}};
}

auto poly_from_autocorr_factors(const std::vector<Vec2>& vrs) -> std::vector<double> {
    if (vrs.empty()) {
        return {1.0};
    }
    // Each factor x^2 - r*x - q contributes 2 roots. For palindromic/autocorrelation
    // polynomials, the reciprocal of each root is also a root. Collect all roots
    // and their reciprocals, then reconstruct with Leja ordering.
    std::vector<std::complex<double>> all_roots;
    all_roots.reserve(4 * vrs.size());
    for (const auto& vr : vrs) {
        auto [r1, r2] = roots_from_quadratic(vr);
        all_roots.push_back(r1);
        all_roots.push_back(r2);
        all_roots.push_back(1.0 / r1);
        all_roots.push_back(1.0 / r2);
    }
    return poly_from_roots(all_roots);
}

void extract_autocorr(Vec2& vr) {
    const auto& r = vr.x();
    const auto& q = vr.y();
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
