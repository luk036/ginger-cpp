#include <algorithm>
#include <cmath>    // for acos, cos, sin
#include <complex>  // for complex, operator*, operator+
#include <ginger/aberth.hpp>
#include <ginger/config.hpp>
#include <lds/lds.hpp>
#include <utility>  // for pair
#include <vector>   // for vector, vector<>::reference, __v...

#include "execution_policy.hpp"  // for ginger::detail::aberth_step, aberth_autocorr_step

using std::vector;
using Complex = std::complex<double>;

// static const auto TWO_PI = 2.0 * std::acos(-1.0);

/**
 * @brief Initial guess for the Aberth-Ehrlich method
 *
 * The `initial_aberth` function calculates the initial values for the Aberth-Ehrlich method for
 * finding the roots of a polynomial.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of a
 * polynomial. Each element of the vector corresponds to a term in the polynomial, starting from the
 * highest degree term and ending with the constant term.
 *
 * @return The function `initial_aberth` returns a vector of Complex numbers representing the
 * initial guesses for the roots of the polynomial.
 */
auto initial_aberth(const vector<double>& coeffs) -> vector<Complex> {
    const auto degree = coeffs.size() - 1;
    const auto center = -coeffs[1] / (static_cast<double>(degree) * coeffs[0]);
    const auto p_center = horner_eval_f(coeffs, center);
    const auto radius = std::pow(std::fabs(p_center), 1.0 / static_cast<double>(degree));
    auto z0s = vector<Complex>{};
    z0s.reserve(degree);
    // lds::Circle<2> c_gen{};
    for (auto i = 0U; i != degree; ++i) {
        // auto res = c_gen.pop();
        auto z0 = center
                  + radius
                        * Complex{
                            circle2_table_y(i),
                            circle2_table_x(
                                i)};  // note! swap x and y to get correct distribution for autocorr
        z0s.emplace_back(z0);
    }
    return z0s;
}

auto aberth(const vector<double>& coeffs, vector<Complex>& zs,
            const ginger::Options& options = ginger::Options()) -> std::pair<unsigned int, bool> {
    ginger::detail::aberth_step step{coeffs};
    return ginger::detail::sequential_policy::run(zs, options, step);
}

/**
 * @brief Initial guess for the Aberth-Ehrlich method (specifically for auto-correlation functions)
 *
 * The `initial_aberth_autocorr` function calculates the initial values for the Aberth-Ehrlich
 * method for finding the roots of a polynomial.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of a
 * polynomial. Each element of the vector corresponds to a term in the polynomial, starting from the
 * highest degree term and ending with the constant term.
 *
 * @return The function `initial_aberth_autocorr` returns a vector of Complex numbers representing
 * the initial guesses for the roots of the polynomial.
 */
auto initial_aberth_autocorr(const vector<double>& coeffs) -> vector<Complex> {
    const auto degree = coeffs.size() - 1;  // assume even
    const auto center = -coeffs[1] / (static_cast<double>(degree) * coeffs[0]);
    const auto poly_c = horner_eval_f(coeffs, center);
    auto radius = std::pow(std::fabs(poly_c), 1.0 / static_cast<double>(degree));
    if (std::abs(radius) > 1.0) {
        radius = 1.0 / radius;
    }
    auto z0s = vector<Complex>{};
    z0s.reserve(degree / 2);
    for (auto i = 0U; i != degree / 2; ++i) {
        auto z0 = center
                  + radius
                        * Complex{circle2_table_y(i),
                                  circle2_table_x(i)};  // note! swap x and y to get correct
                                                        // distribution for autocorr
        z0s.emplace_back(z0);
    }
    return z0s;
}

auto aberth_autocorr(const vector<double>& coeffs, vector<Complex>& zs,
                     const ginger::Options& options = ginger::Options())
    -> std::pair<unsigned int, bool> {
    ginger::detail::aberth_autocorr_step step{coeffs};
    return ginger::detail::sequential_policy::run(zs, options, step);
}

auto leja_order(const vector<Complex>& points) -> vector<Complex> {
    if (points.empty()) {
        return {};
    }
    // Greedy Leja ordering. Maintaining the running minimum distance to the
    // selected set turns the naive O(n^3) rescan into O(n^2).
    auto remaining = points;
    std::sort(remaining.begin(), remaining.end(),
              [](const Complex& a, const Complex& b) { return std::abs(a) < std::abs(b); });
    vector<Complex> result;
    result.reserve(remaining.size());
    result.push_back(remaining.front());
    remaining.erase(remaining.begin());
    vector<double> dists;
    dists.reserve(remaining.size());
    for (const auto& p : remaining) {
        dists.push_back(std::abs(p - result.front()));
    }
    while (!remaining.empty()) {
        auto best_idx = size_t{0};
        auto best_dist = dists.front();
        for (auto i = size_t{1}; i < dists.size(); ++i) {
            if (dists[i] > best_dist) {
                best_dist = dists[i];
                best_idx = i;
            }
        }
        const auto chosen = remaining[best_idx];
        remaining.erase(remaining.begin() + static_cast<ptrdiff_t>(best_idx));
        dists.erase(dists.begin() + static_cast<ptrdiff_t>(best_idx));
        result.push_back(chosen);
        for (auto i = size_t{0}; i < remaining.size(); ++i) {
            const auto dist = std::abs(remaining[i] - chosen);
            if (dist < dists[i]) {
                dists[i] = dist;
            }
        }
    }
    return result;
}

auto poly_from_autocorr_roots(const vector<Complex>& zs) -> vector<double> {
    if (zs.empty()) {
        return {1.0};
    }
    // Add reciprocals to account for the palindromic root-pair structure
    vector<Complex> all_roots;
    all_roots.reserve(2 * zs.size());
    for (const auto& z : zs) {
        all_roots.push_back(z);
        all_roots.push_back(1.0 / z);
    }
    return poly_from_roots(all_roots);
}

auto poly_from_roots(const vector<Complex>& zs) -> vector<double> {
    auto ordered = leja_order(zs);
    vector<Complex> coeffs{1.0};
    for (const auto& z : ordered) {
        auto prev = coeffs[0];
        for (auto i = 1U; i < coeffs.size(); ++i) {
            auto old = coeffs[i];
            coeffs[i] = coeffs[i] - z * prev;
            prev = old;
        }
        coeffs.push_back(-z * prev);
    }
    vector<double> result;
    result.reserve(coeffs.size());
    for (const auto& c : coeffs) {
        result.push_back(c.real());
    }
    return result;
}
