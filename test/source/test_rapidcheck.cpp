#include <doctest/doctest.h>

#include <complex>
#include <vector>

#include "ginger/aberth.hpp"
#include "ginger/config.hpp"
#include "ginger/rootfinding.hpp"

#ifdef RAPIDCHECK_H
#    include <rapidcheck.h>

// Simple helper functions to avoid complex lambdas that cause MSVC ICE
void test_initial_guess_count() {
    auto degree = static_cast<size_t>(*rc::gen::inRange(1, 10));
    std::vector<double> coeffs(degree + 1);
    coeffs[0] = *rc::gen::nonZero<double>();
    for (size_t i = 1; i <= degree; ++i) {
        coeffs[i] = *rc::gen::arbitrary<double>();
    }

    auto initial = initial_aberth(coeffs);

    RC_ASSERT(initial.size() == degree);
}

void test_options_defaults() {
    Options options;
    RC_ASSERT(options.max_iters > static_cast<unsigned int>(0));
    RC_ASSERT(options.tolerance > 0.0);
}

void test_bairstow_initial_guess() {
    auto degree = static_cast<size_t>(*rc::gen::inRange(2, 7) * 2);

    std::vector<double> coeffs(degree + 1);
    coeffs[0] = *rc::gen::nonZero<double>();
    for (size_t i = 1; i <= degree; ++i) {
        coeffs[i] = *rc::gen::arbitrary<double>();
    }

    auto initial = initial_guess(coeffs);

    RC_ASSERT(initial.size() == degree / 2);
}

void test_roots_of_unity() {
    auto n = static_cast<size_t>(*rc::gen::inRange(2, 10));

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

TEST_CASE("Property-based test: Initial guess produces correct number of roots") {
    rc::check("initial_aberth produces degree number of initial guesses", test_initial_guess_count);
}

TEST_CASE("Property-based test: Options default values are reasonable") {
    rc::check("Options can be used with default values", test_options_defaults);
}

TEST_CASE("Property-based test: Bairstow method produces quadratic factors") {
    rc::check("Bairstow initial guess produces correct number of quadratic factors",
              test_bairstow_initial_guess);
}

TEST_CASE("Property-based test: Monic polynomial (x^n - 1) roots are nth roots of unity") {
    rc::check("Roots of x^n - 1 are nth roots of unity", test_roots_of_unity);
}

#endif