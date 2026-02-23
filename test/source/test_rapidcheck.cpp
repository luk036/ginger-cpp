#include <complex>
#include <vector>

#include <doctest/doctest.h>
#include "ginger/aberth.hpp"
#include "ginger/config.hpp"
#include "ginger/rootfinding.hpp"

#ifdef RAPIDCHECK_H
#    include <rapidcheck.h>

// Simple helper functions to avoid complex lambdas that cause MSVC ICE
void test_linear_polynomial() {
    auto a = *rc::gen::nonZero<double>();
    auto b = *rc::gen::arbitrary<double>();
    std::vector<double> coeffs = {a, b};
    
    auto initial = initial_aberth(coeffs);
    Options options;
    std::pair<unsigned int, bool> result = aberth(coeffs, initial, options);
    bool converged = result.second;
    
    RC_ASSERT(converged);
    RC_ASSERT(initial.size() == static_cast<size_t>(1));
}

void test_quadratic_polynomial() {
    auto a = *rc::gen::nonZero<double>();
    auto b = *rc::gen::arbitrary<double>();
    auto c = *rc::gen::arbitrary<double>();
    std::vector<double> coeffs = {a, b, c};
    
    auto initial = initial_aberth(coeffs);
    Options options;
    std::pair<unsigned int, bool> result = aberth(coeffs, initial, options);
    bool converged = result.second;
    
    RC_ASSERT(converged);
    RC_ASSERT(initial.size() == static_cast<size_t>(2));
}

void test_horner_consistency() {
    auto degree = *rc::gen::inRange(2, 10);
    std::vector<double> coeffs;
    coeffs.reserve(degree + 1);
    for (size_t i = 0; i <= degree; ++i) {
        coeffs.push_back(*rc::gen::arbitrary<double>());
    }
    while (coeffs[0] == 0.0) {
        coeffs[0] = *rc::gen::nonZero<double>();
    }
    
    auto x = *rc::gen::arbitrary<double>();
    double horner_value = horner_eval(coeffs, degree, x);
    
    double direct_value = 0.0;
    for (size_t i = 0; i <= degree; ++i) {
        direct_value += coeffs[i] * std::pow(x, static_cast<int>(degree - i));
    }
    
    RC_ASSERT(std::abs(horner_value - direct_value) < 1e-10);
}

void test_initial_guess_count() {
    auto degree = *rc::gen::inRange(1, 10);
    std::vector<double> coeffs(degree + 1);
    coeffs[0] = *rc::gen::nonZero<double>();
    for (size_t i = 1; i <= degree; ++i) {
        coeffs[i] = *rc::gen::arbitrary<double>();
    }
    
    auto initial = initial_aberth(coeffs);
    
    RC_ASSERT(initial.size() == static_cast<size_t>(degree));
}

void test_options_defaults() {
    Options options;
    RC_ASSERT(options.max_iters > static_cast<unsigned int>(0));
    RC_ASSERT(options.tolerance > 0.0);
}

void test_bairstow_initial_guess() {
    auto degree = *rc::gen::inRange(2, 7) * 2;
    
    std::vector<double> coeffs(degree + 1);
    coeffs[0] = *rc::gen::nonZero<double>();
    for (size_t i = 1; i <= degree; ++i) {
        coeffs[i] = *rc::gen::arbitrary<double>();
    }
    
    auto initial = initial_guess(coeffs);
    
    RC_ASSERT(initial.size() == degree / 2);
}

void test_roots_of_unity() {
    auto n = *rc::gen::inRange(2, 10);
    
    std::vector<double> coeffs(n + 1);
    coeffs[0] = 1.0;
    for (size_t i = 1; i < n; ++i) {
        coeffs[i] = 0.0;
    }
    coeffs[n] = -1.0;
    
    auto initial = initial_aberth(coeffs);
    Options options;
    std::pair<unsigned int, bool> result = aberth(coeffs, initial, options);
    bool converged = result.second;
    
    if (converged) {
        for (const auto& root : initial) {
            std::complex<double> z = std::pow(root, static_cast<int>(n));
            RC_ASSERT(std::abs(z - 1.0) < options.tolerance);
        }
    }
}

void test_autocorr_conjugate_pairs() {
    auto degree = *rc::gen::inRange(2, 6);
    
    std::vector<double> coeffs(degree + 1);
    coeffs[0] = 1.0;
    for (size_t i = 1; i <= degree; ++i) {
        coeffs[i] = *rc::gen::arbitrary<double>();
    }
    coeffs[degree] = 1.0;
    
    auto initial = initial_aberth_autocorr(coeffs);
    Options options;
    std::pair<unsigned int, bool> result = aberth_autocorr(coeffs, initial, options);
    bool converged = result.second;
    
    if (converged) {
        for (const auto& root : initial) {
            if (std::abs(root.imag()) > options.tolerance) {
                bool found_conjugate = false;
                for (const auto& other : initial) {
                    if (std::abs(root - std::conj(other)) < options.tolerance) {
                        found_conjugate = true;
                        break;
                    }
                }
                RC_ASSERT(found_conjugate);
            }
        }
    }
}

void test_polynomial_at_roots() {
    auto degree = *rc::gen::inRange(2, 5);
    std::vector<double> coeffs;
    coeffs.reserve(degree + 1);
    for (size_t i = 0; i <= degree; ++i) {
        coeffs.push_back(*rc::gen::arbitrary<double>());
    }
    while (coeffs[0] == 0.0) {
        coeffs[0] = *rc::gen::nonZero<double>();
    }
    
    auto initial = initial_aberth(coeffs);
    Options options;
    std::pair<unsigned int, bool> result = aberth(coeffs, initial, options);
    bool converged = result.second;
    
    if (converged) {
        for (const auto& root : initial) {
            std::complex<double> z(root);
            std::complex<double> value_complex(0.0);
            for (size_t i = 0; i < coeffs.size(); ++i) {
                value_complex = value_complex * z + std::complex<double>(coeffs[i]);
            }
            RC_ASSERT(std::abs(value_complex) < options.tolerance);
        }
    }
}

TEST_CASE("Property-based test: Aberth method finds correct roots for linear polynomial") {
    rc::check("aberth finds root of ax + b", test_linear_polynomial);
}

TEST_CASE("Property-based test: Aberth method finds correct roots for quadratic polynomial") {
    rc::check("aberth finds roots of ax^2 + bx + c", test_quadratic_polynomial);
}

TEST_CASE("Property-based test: Horner evaluation consistency") {
    rc::check("horner_eval is consistent with direct evaluation", test_horner_consistency);
}

TEST_CASE("Property-based test: Initial guess produces correct number of roots") {
    rc::check("initial_aberth produces degree number of initial guesses", test_initial_guess_count);
}

TEST_CASE("Property-based test: Options default values are reasonable") {
    rc::check("Options can be used with default values", test_options_defaults);
}

TEST_CASE("Property-based test: Bairstow method produces quadratic factors") {
    rc::check("Bairstow initial guess produces correct number of quadratic factors", test_bairstow_initial_guess);
}

TEST_CASE("Property-based test: Monic polynomial (x^n - 1) roots are nth roots of unity") {
    rc::check("Roots of x^n - 1 are nth roots of unity", test_roots_of_unity);
}

TEST_CASE("Property-based test: Aberth autocorr finds conjugate pairs") {
    rc::check("aberth_autocorr finds conjugate pairs for autocorrelation polynomials", test_autocorr_conjugate_pairs);
}

TEST_CASE("Property-based test: Polynomial evaluation at roots is near zero") {
    rc::check("Evaluating polynomial at found roots gives near-zero values", test_polynomial_at_roots);
}

#endif