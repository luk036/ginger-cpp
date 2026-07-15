#include <algorithm>
#include <cmath>    // for acos, cos, sin
#include <complex>  // for complex, operator*, operator+
#include <future>   // for future
#include <ginger/config.hpp>
#include <ginger/robin.hpp>        // for Robin
#include <ginger/thread_pool.hpp>  // for thread_pool
#include <lds/lds.hpp>
#include <limits>   // for numeric_limits
#include <utility>  // for pair
#include <vector>   // for vector, vector<>::reference, __v...

using std::vector;
using Complex = std::complex<double>;

// static const auto TWO_PI = 2.0 * std::acos(-1.0);

/// @brief Helper to generate a constexpr table of VdCorput<Base> values
/// @tparam N Number of values to generate
/// @tparam Base Base of the van der Corput sequence
/// @return std::array<double, N> with precomputed sequence values
template <unsigned long N, unsigned long Base = 2> constexpr auto make_vdc_table()
    -> std::array<double, N> {
    std::array<double, N> table{};
    lds::VdCorput<Base> gen;
    for (unsigned long i = 0; i < N; ++i) {
        table[i] = gen.pop();
    }
    return table;
}

/// @brief Size of the precomputed VdCorput base-2 table
constexpr const auto VDC_TABLE_SIZE = 1000UL;

/// @brief Precomputed table of VdCorput sequence values (base 2)
/// @details Generated at compile-time using VdCorput<2>
constexpr std::array<double, VDC_TABLE_SIZE> VDC_TABLE_2 = make_vdc_table<VDC_TABLE_SIZE, 2>();

/// @brief Access the precomputed VdCorput base-2 table
/// @param index Index into the table
/// @return The VDC value at the given index
double vdc2_table(unsigned long index) { return VDC_TABLE_2[index]; }

/// @brief Precomputed table of 1000 Circle<2> points
/// @details Generated using precomputed VDC_TABLE_2 mapped to unit circle.
///          Not constexpr: std::cos/std::sin lack portable constexpr support in C++20.
static const auto CIRCLE_TABLE_2 = []() {
    std::array<std::array<double, 2>, VDC_TABLE_SIZE> table{};
    for (unsigned long i = 0; i < VDC_TABLE_SIZE; ++i) {
        auto theta = VDC_TABLE_2[i] * lds::TWO_PI;
        table[i] = {std::cos(theta), std::sin(theta)};
    }
    return table;
}();

/// @brief Access the precomputed Circle base-2 table x-coordinate
/// @param index Index into the table
/// @return The circle point x-coordinate at the given index
constexpr double circle2_table_x(unsigned long index) { return CIRCLE_TABLE_2[index][0]; }

/// @brief Access the precomputed Circle base-2 table y-coordinate
/// @param index Index into the table
/// @return The circle point y-coordinate at the given index
constexpr double circle2_table_y(unsigned long index) { return CIRCLE_TABLE_2[index][1]; }

/// @brief Precomputed table of cos(pi * vdc2_table[i]) values
/// @details Used by initial_guess and initial_autocorr to avoid computing cos on the fly.
static const auto COS_PI_VDC2_TABLE = []() {
    std::array<double, VDC_TABLE_SIZE> table{};
    for (unsigned long i = 0; i < VDC_TABLE_SIZE; ++i) {
        table[i] = std::cos(lds::TWO_PI / 2.0 * VDC_TABLE_2[i]);
    }
    return table;
}();

/// @brief Access the precomputed cos(pi * vdc2_table[i]) value
/// @param index Index into the table
/// @return The cos(pi * vdc_value) at the given index
double cos_pi_vdc2(unsigned long index) { return COS_PI_VDC2_TABLE[index]; }

/**
 * The function `horner_eval_c` is implementing the Horner's method for
 * evaluating a polynomial at a given point.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of a
 * polynomial. Each element of the vector corresponds to a term in the polynomial, starting from the
 * highest degree term and ending with the constant term.
 * @param[in] z
 * @return Tp
 *
 * @verbatim
 *        coeffs[0]      coeffs[1]      coeffs[2]                 coeffs[n-1]    coeffs[n]
 *     +-----------> + -----------> + -----------> + ... + -----------> + --------->
 *     |             |              |              |                   |            |
 *     |    z       v    z         v    z         v                   v    z       v    z
 *     +---> [x] --> +---> [x] --> +---> [x] --> +---> ... ----> +---> [x] --> +---> [x] --> result
 *           |                    |              |                           |              |
 *           +---------------------              +---------------------------+              |
 *           |                                   |                                        |
 *           +------------------------------------+----------------------------------------+
 *
 * P(x) = coeffs[0]*x^n + coeffs[1]*x^(n-1) + ... + coeffs[n-1]*x + coeffs[n]
 * @endverbatim
 */
inline auto horner_eval_c(const std::vector<double>& coeffs, const std::complex<double>& zval)
    -> std::complex<double> {
    std::complex<double> result(0.0, 0.0);
    for (auto coeff : coeffs) {
        result = result * zval + coeff;
    }
    return result;
}

/**
 * The function `horner_eval_f` is implementing the Horner's method for
 * evaluating a polynomial at a given point.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of a
 * polynomial. Each element of the vector corresponds to a term in the polynomial, starting from the
 * highest degree term and ending with the constant term.
 * @param[in] z
 * @return Tp
 */
inline auto horner_eval_f(const std::vector<double>& coeffs, const double& zval) -> double {
    double result(0.0);
    for (auto coeff : coeffs) {
        result = result * zval + coeff;
    }
    return result;
}

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
 *
 * @verbatim
 *        center
 *          *
 *         /|\
 *        / | \ radius
 *       /  |  \
 *      *   |   *
 *     /    |    \
 *    *----------*----------> real
 *     \    |    /
 *      *   |   *
 *       \  |  /
 *        \ | /
 *         \|/
 *          *
 *         imag
 *
 * Initial points distributed on a circle around center
 * @endverbatim
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

// ST core — no std::future overhead, returns tolerance directly
template <typename F> static auto aberth_st_core(const vector<double>& coeffs, vector<Complex>& zs,
                                                 const Options& options, F& aberth_job_generator)
    -> std::pair<unsigned int, bool> {
    const auto num_roots = zs.size();
    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        auto aberth_job = aberth_job_generator(coeffs, zs);
        for (auto idx = 0U; idx != num_roots; ++idx) {
            tolerance = std::max(tolerance, aberth_job(idx));
        }
        if (tolerance < options.tolerance) {
            return {niter, true};
        }
    }
    return {options.max_iters, false};
}

// MT core — uses futures from thread pool
template <typename F> static auto aberth_mt_core(const vector<double>& coeffs, vector<Complex>& zs,
                                                 const Options& options, ginger::thread_pool& pool,
                                                 F& aberth_job_generator)
    -> std::pair<unsigned int, bool> {
    const auto num_roots = zs.size();
    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        auto aberth_job = aberth_job_generator(coeffs, zs);
        vector<std::future<double>> results;
        results.reserve(num_roots);
        for (auto idx = 0U; idx != num_roots; ++idx) {
            results.emplace_back(pool.enqueue([&, idx]() { return aberth_job(idx); }));
        }
        for (auto& result : results) {
            tolerance = std::max(tolerance, result.get());
        }
        if (tolerance < options.tolerance) {
            return {niter, true};
        }
    }
    return {options.max_iters, false};
}

auto aberth(const vector<double>& coeffs, vector<Complex>& zs, const Options& options = Options())
    -> std::pair<unsigned int, bool> {
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto idx = 0U; idx != degree; ++idx) {
        coeffs1[idx] = static_cast<double>(degree - idx) * coeffs[idx];
    }
    const auto num_zs = zs.size();

    auto aberth_job_generator = [&](const vector<double>&, vector<Complex>& zs_ref) {
        return [&, num_zs](size_t idx) -> double {
            const auto zi = zs_ref[idx];
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            for (auto jdx = 0U; jdx < num_zs; ++jdx) {
                if (jdx == idx) continue;
                P1 -= P / (zi - zs_ref[jdx]);
            }
            zs_ref[idx] -= P / P1;
            return tol_i;
        };
    };

    return aberth_st_core(coeffs, zs, options, aberth_job_generator);
}

auto aberth_mt(const vector<double>& coeffs, vector<Complex>& zs,
               const Options& options = Options()) -> std::pair<unsigned int, bool> {
    auto& pool = ginger::get_thread_pool();
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto idx = 0U; idx != degree; ++idx) {
        coeffs1[idx] = static_cast<double>(degree - idx) * coeffs[idx];
    }
    const auto num_zs = zs.size();

    auto aberth_job_generator = [&](const vector<double>&, vector<Complex>& zs_ref) {
        return [&, num_zs](size_t idx) -> double {
            const auto zi = zs_ref[idx];
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            for (auto jdx = 0U; jdx < num_zs; ++jdx) {
                if (jdx == idx) continue;
                P1 -= P / (zi - zs_ref[jdx]);
            }
            zs_ref[idx] -= P / P1;
            return tol_i;
        };
    };

    return aberth_mt_core(coeffs, zs, options, pool, aberth_job_generator);
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
 *
 * @verbatim
 * For auto-correlation functions:
 *
 *        center
 *          *
 *         /|\
 *        / | \ radius (limited to 1/radius if > 1)
 *       /  |  \
 *      *   |   *
 *     /    |    \
 *    *----------*----------> real
 *     \    |    /
 *      *   |   *
 *       \  |  /
 *        \ | /
 *         \|/
 *          *
 *         imag
 *
 * Initial points distributed on a circle around center, with radius adjustment
 * for auto-correlation specific properties
 * @endverbatim
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

// ST core — no std::future overhead, returns tolerance directly
template <typename F>
static auto aberth_autocorr_st_core(const vector<double>& coeffs, vector<Complex>& zs,
                                    const Options& options, F& aberth_job_generator)
    -> std::pair<unsigned int, bool> {
    const auto num_roots = zs.size();
    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        auto aberth_job = aberth_job_generator(coeffs, zs);
        for (auto idx = 0U; idx != num_roots; ++idx) {
            tolerance = std::max(tolerance, aberth_job(idx));
        }
        if (tolerance < options.tolerance) {
            return {niter, true};
        }
    }
    return {options.max_iters, false};
}

// MT core — uses futures from thread pool
template <typename F>
static auto aberth_autocorr_mt_core(const vector<double>& coeffs, vector<Complex>& zs,
                                    const Options& options, ginger::thread_pool& pool,
                                    F& aberth_job_generator) -> std::pair<unsigned int, bool> {
    const auto num_roots = zs.size();
    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        auto aberth_job = aberth_job_generator(coeffs, zs);
        vector<std::future<double>> results;
        results.reserve(num_roots);
        for (auto idx = 0U; idx != num_roots; ++idx) {
            results.emplace_back(pool.enqueue([&, idx]() { return aberth_job(idx); }));
        }
        for (auto& result : results) {
            tolerance = std::max(tolerance, result.get());
        }
        if (tolerance < options.tolerance) {
            return {niter, true};
        }
    }
    return {options.max_iters, false};
}

auto aberth_autocorr(const vector<double>& coeffs, vector<Complex>& zs,
                     const Options& options = Options()) -> std::pair<unsigned int, bool> {
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto idx = 0U; idx != degree; ++idx) {
        coeffs1[idx] = static_cast<double>(degree - idx) * coeffs[idx];
    }
    const auto num_zs = zs.size();

    auto aberth_job_generator = [&](const vector<double>&, vector<Complex>& zs_ref) {
        return [&, num_zs](size_t idx) -> double {
            const auto zi = zs_ref[idx];
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            for (auto jdx = 0U; jdx < num_zs; ++jdx) {
                if (jdx == idx) continue;
                P1 -= P / (zi - zs_ref[jdx]);
                P1 -= P / (zi - 1.0 / zs_ref[jdx]);
            }
            zs_ref[idx] -= P / P1;
            return tol_i;
        };
    };

    return aberth_autocorr_st_core(coeffs, zs, options, aberth_job_generator);
}

auto aberth_autocorr_mt(const vector<double>& coeffs, vector<Complex>& zs,
                        const Options& options = Options()) -> std::pair<unsigned int, bool> {
    auto& pool = ginger::get_thread_pool();
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto idx = 0U; idx != degree; ++idx) {
        coeffs1[idx] = static_cast<double>(degree - idx) * coeffs[idx];
    }
    const auto num_zs = zs.size();

    auto aberth_job_generator = [&](const vector<double>&, vector<Complex>& zs_ref) {
        return [&, num_zs](size_t idx) -> double {
            const auto zi = zs_ref[idx];
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            for (auto jdx = 0U; jdx < num_zs; ++jdx) {
                if (jdx == idx) continue;
                P1 -= P / (zi - zs_ref[jdx]);
                P1 -= P / (zi - 1.0 / zs_ref[jdx]);
            }
            zs_ref[idx] -= P / P1;
            return tol_i;
        };
    };

    return aberth_autocorr_mt_core(coeffs, zs, options, pool, aberth_job_generator);
}

auto leja_order(const vector<Complex>& points) -> vector<Complex> {
    if (points.empty()) {
        return {};
    }
    // Greedy Leja ordering (O(n^2)):
    // 1. Start with the smallest-magnitude point
    // 2. Each subsequent point maximizes the MIN distance
    //    to all already-selected points
    auto sorted = points;
    std::sort(sorted.begin(), sorted.end(),
              [](const Complex& a, const Complex& b) { return std::abs(a) < std::abs(b); });
    vector<Complex> result;
    result.reserve(sorted.size());
    result.push_back(sorted.front());
    sorted.erase(sorted.begin());
    while (!sorted.empty()) {
        auto best_idx = size_t{0};
        auto best_dist = -1.0;
        for (auto i = size_t{0}; i < sorted.size(); ++i) {
            auto min_dist = std::numeric_limits<double>::max();
            for (const auto& p : result) {
                min_dist = std::min(min_dist, std::abs(sorted[i] - p));
            }
            if (min_dist > best_dist) {
                best_dist = min_dist;
                best_idx = i;
            }
        }
        result.push_back(sorted[best_idx]);
        sorted.erase(sorted.begin() + static_cast<ptrdiff_t>(best_idx));
    }
    return result;
}

auto poly_from_roots(const vector<Complex>& zs) -> vector<double>;

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
