// Order-independence tests for the multi-threaded (Jacobi) Bairstow variant.
// Mirrors the Rust ginger-rs suite: the Jacobi variant must be order-independent
// (same iteration count and same root SET under permuted initial guesses).
//
// NOTE: pbairstow_even_mt falls back to Gauss-Seidel when vrs.size() <= 4,
// so these tests use the degree-12 palindromic polynomial (6 factors) which
// exercises the true multi-threaded Jacobi snapshot path.
#include <doctest/doctest.h>

#include <algorithm>                  // for sort
#include <complex>                    // for complex
#include <ginger/config.hpp>          // for Options
#include <ginger/rootfinding_mt.hpp>  // for initial_guess, pbairstow_even_mt
#include <vector>                     // for vector

using namespace ginger;

// The roots_from_quadratic helper is static in source/rootfinding.cpp, so it
// is duplicated here (same pattern as examples/order_experiment.rs).
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

static auto sorted_roots(const std::vector<Vec2>& vrs) -> std::vector<std::pair<double, double>> {
    std::vector<std::pair<double, double>> roots;
    for (const auto& vr : vrs) {
        auto [a, b] = roots_from_quadratic(vr);
        roots.emplace_back(a.real(), a.imag());
        roots.emplace_back(b.real(), b.imag());
    }
    std::sort(roots.begin(), roots.end());
    return roots;
}

static auto max_root_set_diff(const std::vector<std::pair<double, double>>& a,
                              const std::vector<std::pair<double, double>>& b) -> double {
    auto worst = 0.0;
    for (auto i = 0U; i < a.size(); ++i) {
        worst = std::max(worst, std::abs(a[i].first - b[i].first));
        worst = std::max(worst, std::abs(a[i].second - b[i].second));
    }
    return worst;
}

TEST_CASE("test jacobi mt order independent") {
    // Palindromic degree-12: 6 factors, exercises the true Jacobi path.
    const auto h
        = std::vector<double>{1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 3.0, 0.0, 2.0, 0.0, 1.0};
    auto options = Options();
    options.tolerance = 1e-12;

    const auto base = initial_guess(h);
    REQUIRE_EQ(base.size(), 6U);

    std::vector<std::vector<Vec2>> perms{
        base,
        {base[5], base[4], base[3], base[2], base[1], base[0]},
        {base[1], base[4], base[0], base[5], base[2], base[3]},
        {base[3], base[1], base[5], base[0], base[4], base[2]},
    };

    std::vector<unsigned int> niters;
    std::vector<std::vector<std::pair<double, double>>> rootsets;
    for (auto& perm : perms) {
        auto result = pbairstow_even_mt(h, perm, options);
        REQUIRE(result.second);
        niters.push_back(result.first);
        rootsets.push_back(sorted_roots(perm));
    }

    // Same iteration count for every permutation.
    for (const auto n : niters) {
        CHECK_EQ(n, niters[0]);
    }
    // Same converged root SET (up to floating-point noise).
    for (auto k = 1U; k < rootsets.size(); ++k) {
        CHECK_LT(max_root_set_diff(rootsets[0], rootsets[k]), 1e-9);
    }
}
