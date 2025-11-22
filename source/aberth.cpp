#include <ginger/ThreadPool.h>  // for ThreadPool

#include <cmath>    // for acos, cos, sin
#include <complex>  // for complex, operator*, operator+
#include <future>   // for future
#include <ginger/config.hpp>
#include <ginger/robin.hpp>  // for Robin
#include <ldsgen/lds.hpp>
#include <utility>  // for pair
#include <vector>   // for vector, vector<>::reference, __v...

using std::vector;
using Complex = std::complex<double>;

// static const auto TWO_PI = 2.0 * std::acos(-1.0);

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
 * ```svgbob
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
 * ```
 */
inline auto horner_eval_c(const std::vector<double> &coeffs, const std::complex<double> &zval)
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
inline auto horner_eval_f(const std::vector<double> &coeffs, const double &zval) -> double {
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
 * ```svgbob
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
 * ```
 */
auto initial_aberth(const vector<double> &coeffs) -> vector<Complex> {
    const auto degree = coeffs.size() - 1;
    const auto center = -coeffs[1] / (double(degree) * coeffs[0]);
    const auto p_center = horner_eval_f(coeffs, center);
    const auto radius = std::pow(std::fabs(p_center), 1.0 / double(degree));
    auto z0s = vector<Complex>{};
    auto c_gen = ldsgen::Circle(2);
    for (auto i = 0U; i != degree; ++i) {
        auto res = c_gen.pop();
        auto z0 = center + radius * Complex{res[1], res[0]};
        z0s.emplace_back(z0);
    }
    return z0s;
}

template <typename F>
auto aberth_core(const vector<double> &coeffs, vector<Complex> &zs, const Options &options, F &&aberth_job_generator)
    -> std::pair<unsigned int, bool> {
    const auto m = zs.size();
    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        auto aberth_job = aberth_job_generator(coeffs, zs);
        vector<std::future<double>> results;
        for (auto i = 0U; i != m; ++i) {
            results.emplace_back(aberth_job(i));
        }
        for (auto &result : results) {
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

auto aberth(const vector<double> &coeffs, vector<Complex> &zs, const Options &options = Options())
    -> std::pair<unsigned int, bool> {
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto i = 0U; i != degree; ++i) {
        coeffs1[i] = double(degree - i) * coeffs[i];
    }
    const auto rr = fun::Robin<size_t>(zs.size());

    auto aberth_job_generator = [&](const vector<double>&, vector<Complex>& zs_ref) {
        return [&](size_t i) {
            const auto zi = zs_ref[i];
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            for (auto j : rr.exclude(i)) {
                P1 -= P / (zi - zs_ref[j]);
            }
            zs_ref[i] -= P / P1;
            return std::async(std::launch::deferred, [tol_i](){ return tol_i; });
        };
    };

    return aberth_core(coeffs, zs, options, aberth_job_generator);
}

auto aberth_mt(const vector<double> &coeffs, vector<Complex> &zs,
               const Options &options = Options()) -> std::pair<unsigned int, bool> {
    ThreadPool pool(std::thread::hardware_concurrency());
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto i = 0U; i != degree; ++i) {
        coeffs1[i] = double(degree - i) * coeffs[i];
    }
    const auto rr = fun::Robin<size_t>(zs.size());

    auto aberth_job_generator = [&](const vector<double>&, vector<Complex>& zs_ref) {
        return [&](size_t i) {
            return pool.enqueue([&, i]() {
                const auto zi = zs_ref[i];
                const auto P = horner_eval_c(coeffs, zi);
                const auto tol_i = std::abs(P);
                auto P1 = horner_eval_c(coeffs1, zi);
                for (auto j : rr.exclude(i)) {
                    P1 -= P / (zi - zs_ref[j]);
                }
                zs_ref[i] -= P / P1;
                return tol_i;
            });
        };
    };

    return aberth_core(coeffs, zs, options, aberth_job_generator);
}

/**
 * @brief Initial guess for the Aberth-Ehrlich method (specifically for auto-correlation functions)
 *
 * The `initial_aberth_autocorr` function calculates the initial values for the Aberth-Ehrlich method for
 * finding the roots of a polynomial.
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of a
 * polynomial. Each element of the vector corresponds to a term in the polynomial, starting from the
 * highest degree term and ending with the constant term.
 *
 * @return The function `initial_aberth_autocorr` returns a vector of Complex numbers representing the
 * initial guesses for the roots of the polynomial.
 *
 * ```svgbob
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
 * ```
 */
auto initial_aberth_autocorr(const vector<double> &coeffs) -> vector<Complex> {
    const auto degree = coeffs.size() - 1;  // assume even
    const auto center = -coeffs[1] / (double(degree) * coeffs[0]);
    const auto poly_c = horner_eval_f(coeffs, center);
    auto radius = std::pow(std::fabs(poly_c), 1.0 / double(degree));
    if (std::abs(radius) > 1.0) {
        radius = 1.0 / radius;
    }
    auto z0s = vector<Complex>{};
    auto c_gen = ldsgen::Circle(2);
    for (auto i = 0U; i != degree / 2; ++i) {
        auto res = c_gen.pop();
        auto z0 = center + radius * Complex{res[1], res[0]};
        z0s.emplace_back(z0);
    }
    return z0s;
}

template <typename F>
auto aberth_autocorr_core(const vector<double> &coeffs, vector<Complex> &zs, const Options &options, F &&aberth_job_generator)
    -> std::pair<unsigned int, bool> {
    const auto m = zs.size();
    for (auto niter = 0U; niter != options.max_iters; ++niter) {
        auto tolerance = 0.0;
        auto aberth_job = aberth_job_generator(coeffs, zs);
        vector<std::future<double>> results;
        for (auto i = 0U; i != m; ++i) {
            results.emplace_back(aberth_job(i));
        }
        for (auto &result : results) {
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

auto aberth_autocorr(const vector<double> &coeffs, vector<Complex> &zs,
                     const Options &options = Options()) -> std::pair<unsigned int, bool> {
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto i = 0U; i != degree; ++i) {
        coeffs1[i] = double(degree - i) * coeffs[i];
    }
    const auto rr = fun::Robin<size_t>(zs.size());

    auto aberth_job_generator = [&](const vector<double>&, vector<Complex>& zs_ref) {
        return [&](size_t i) {
            const auto zi = zs_ref[i];
            const auto P = horner_eval_c(coeffs, zi);
            const auto tol_i = std::abs(P);
            auto P1 = horner_eval_c(coeffs1, zi);
            for (auto j : rr.exclude(i)) {
                P1 -= P / (zi - zs_ref[j]);
                P1 -= P / (zi - 1.0 / zs_ref[j]);
            }
            zs_ref[i] -= P / P1;
            return std::async(std::launch::deferred, [tol_i](){ return tol_i; });
        };
    };

    return aberth_autocorr_core(coeffs, zs, options, aberth_job_generator);
}

auto aberth_autocorr_mt(const vector<double> &coeffs, vector<Complex> &zs,
                        const Options &options = Options()) -> std::pair<unsigned int, bool> {
    ThreadPool pool(std::thread::hardware_concurrency());
    const auto degree = coeffs.size() - 1;
    auto coeffs1 = vector<double>(degree);
    for (auto i = 0U; i != degree; ++i) {
        coeffs1[i] = double(degree - i) * coeffs[i];
    }
    const auto rr = fun::Robin<size_t>(zs.size());

    auto aberth_job_generator = [&](const vector<double>&, vector<Complex>& zs_ref) {
        return [&](size_t i) {
            return pool.enqueue([&, i]() {
                const auto zi = zs_ref[i];
                const auto P = horner_eval_c(coeffs, zi);
                const auto tol_i = std::abs(P);
                auto P1 = horner_eval_c(coeffs1, zi);
                for (auto j : rr.exclude(i)) {
                    P1 -= P / (zi - zs_ref[j]);
                    P1 -= P / (zi - 1.0 / zs_ref[j]);
                }
                zs_ref[i] -= P / P1;
                return tol_i;
            });
        };
    };

    return aberth_autocorr_core(coeffs, zs, options, aberth_job_generator);
}