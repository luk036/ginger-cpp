#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <future>
#include <vector>

#include "./include/ginger/ThreadPool.h"  // Include the ThreadPool implementation discussed earlier

static const double TWO_PI = 6.283185307179586476925286766559;  // M_PI * 2

// Assuming Options struct exists and is defined elsewhere
struct Options {
    size_t max_iters;
    double tolerance;
};

auto horner_eval_f(const std::vector<double>& coeffs, double zval) -> double {
    double result = 0.0;
    for (const auto& coeff : coeffs) {
        result = result * zval + coeff;
    }
    return result;
}

auto horner_eval_c(const std::vector<double>& coeffs,
                   const std::complex<double>& zval) -> std::complex<double> {
    std::complex<double> result(0.0, 0.0);
    for (auto coeff : coeffs) {
        result = result * zval + coeff;
    }
    return result;
}

std::vector<std::complex<double>> initial_aberth(const std::vector<double>& coeffs) {
    size_t degree = coeffs.size() - 1;
    double center = -coeffs[1] / (coeffs[0] * degree);
    double poly_c = horner_eval_f(coeffs, center);
    std::complex<double> radius(-poly_c, 0.0);
    radius = std::pow(radius, 1.0 / degree);
    double k = TWO_PI / degree;
    std::vector<std::complex<double>> roots;
    for (size_t idx = 0; idx < degree; ++idx) {
        double theta = k * (0.25 + idx);
        roots.push_back(center + radius * std::complex<double>(cos(theta), sin(theta)));
    }
    return roots;
}

double aberth_job(const std::vector<double>& coeffs, size_t i, std::complex<double>& zi,
                  const std::vector<std::complex<double>>& zsc,
                  const std::vector<double>& coeffs1) {
    std::complex<double> pp = horner_eval_c(coeffs, zi);
    double tol_i = std::abs(pp);  // Using l1_norm equivalent for simplicity
    std::complex<double> pp1 = horner_eval_c(coeffs1, zi);
    for (size_t j = 0; j < zsc.size(); ++j) {
        if (i != j) {
            pp1 -= pp / (zi - zsc[j]);
        }
    }
    zi -= pp / pp1;
    return tol_i;
}

std::pair<int, bool> aberth(const std::vector<double>& coeffs,
                            std::vector<std::complex<double>>& zs, const Options& options) {
    size_t m_zs = zs.size();
    size_t degree = coeffs.size() - 1;
    std::vector<double> coeffs1(degree);
    for (size_t i = 0; i < degree; ++i) {
        coeffs1[i] = coeffs[i] * (degree - i);
    }

    for (size_t niter = 0; niter < options.max_iters; ++niter) {
        double tolerance = 0.0;

        for (size_t i = 0; i < m_zs; ++i) {
            std::complex<double> zi = zs[i];
            double tol_i = aberth_job(coeffs, i, zi, zs, coeffs1);
            if (tol_i > tolerance) tolerance = tol_i;
            zs[i] = zi;
        }
        if (tolerance < options.tolerance) return {niter, true};
    }
    return {options.max_iters, false};
}

// Note: The multi-threaded version using Rayon in Rust would require a different approach in C++,
// potentially using std::thread, OpenMP, or another threading library, which significantly
// increases complexity and is thus omitted here.

// Tests would be implemented with a testing framework like Google Test./ ... Other helper functions
// like horner_eval_f, horner_eval_c, initial_aberth, aberth_job ...

std::pair<int, bool> aberth_mt(const std::vector<double>& coeffs,
                               std::vector<std::complex<double>>& zs, const Options& options) {
    ThreadPool pool(std::thread::hardware_concurrency());  // Initialize thread pool

    size_t m_zs = zs.size();
    size_t degree = coeffs.size() - 1;
    std::vector<double> coeffs1(degree);
    for (size_t i = 0; i < degree; ++i) {
        coeffs1[i] = coeffs[i] * (degree - i);
    }
    std::vector<bool> converged(m_zs, false);

    for (size_t niter = 0; niter < options.max_iters; ++niter) {
        double tolerance = 0.0;

        std::vector<std::future<double>> futures;
        std::vector<std::complex<double>> zsc(zs);  // Copy zs to zsc for thread safety

        for (size_t i = 0; i < zs.size(); ++i) {
            if (!converged[i]) {
                futures.push_back(pool.enqueue(aberth_job, coeffs, i, std::ref(zs[i]),
                                               std::ref(converged[i]), zsc, coeffs1));
            }
        }

        for (auto& future : futures) {
            if (future.valid()) {
                double tol_i = future.get();
                if (tol_i > tolerance) tolerance = tol_i;
            }
        }

        if (tolerance < options.tolerance) return {niter, true};  // Convergence achieved
    }
    return {options.max_iters, false};  // Max iterations reached without convergence
}
