#include <doctest/doctest.h>

#include <chrono>
#include <cmath>
#include <ginger/aberth.hpp>
#include <ginger/autocorr.hpp>
#include <ginger/config.hpp>
#include <ginger/rootfinding.hpp>
#include <random>
#include <vector>

using Vec2 = ginger::Vector2<double>;

/// @brief Generate a random polynomial of given degree with coefficients in [-10, 10]
static auto random_poly(int degree) -> std::vector<double> {
    std::vector<double> h(static_cast<size_t>(degree) + 1);
    std::mt19937_64 rng(
        static_cast<unsigned long>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    for (int i = 0; i <= degree; ++i) {
        h[static_cast<size_t>(i)] = dist(rng);
    }
    return h;
}

/// @brief Generate a random palindromic polynomial (a_i = a_{n-i})
static auto random_palindromic_poly(int degree) -> std::vector<double> {
    REQUIRE(degree % 2 == 0);
    auto h = random_poly(degree / 2);
    // Mirror the first half (excluding the last) to make it palindromic
    for (int i = degree / 2 - 1; i >= 0; --i) {
        h.push_back(h[static_cast<size_t>(i)]);
    }
    return h;
}

// =====================================================================
// Aberth stress tests
// =====================================================================

TEST_CASE("stress test aberth (ST) with high-degree polynomial") {
    const auto h = random_poly(100);

    Options options;
    options.tolerance = 1e-9;
    options.max_iters = 2000;

    auto zs = initial_aberth(h);
    REQUIRE_EQ(zs.size(), h.size() - 1);

    auto result = aberth(h, zs, options);
    auto niter = result.first;
    auto converged = result.second;

    MESSAGE("aberth ST degree-100: niter=", niter, " converged=", converged);
    CHECK_LE(niter, options.max_iters);
    CHECK_GT(niter, 0);
}

TEST_CASE("stress test aberth_mt with high-degree polynomial") {
    const auto h = random_poly(100);

    Options options;
    options.tolerance = 1e-9;
    options.max_iters = 2000;

    auto zs = initial_aberth(h);

    auto result = aberth_mt(h, zs, options);
    auto niter = result.first;
    auto converged = result.second;

    MESSAGE("aberth MT degree-100: niter=", niter, " converged=", converged);
    CHECK_LE(niter, options.max_iters);
    CHECK_GT(niter, 0);
}

TEST_CASE("stress test aberth_autocorr (ST) with high-degree palindromic") {
    const auto h = random_palindromic_poly(100);

    Options options;
    options.tolerance = 1e-9;
    options.max_iters = 2000;

    auto zs = initial_aberth_autocorr(h);
    REQUIRE_EQ(zs.size(), h.size() / 2);

    auto result = aberth_autocorr(h, zs, options);
    auto niter = result.first;
    auto converged = result.second;

    MESSAGE("aberth_autocorr ST degree-100: niter=", niter, " converged=", converged);
    CHECK_LE(niter, options.max_iters);
    CHECK_GT(niter, 0);
}

TEST_CASE("stress test aberth_autocorr (MT) with high-degree palindromic") {
    const auto h = random_palindromic_poly(100);

    Options options;
    options.tolerance = 1e-9;
    options.max_iters = 2000;

    auto zs = initial_aberth_autocorr(h);
    REQUIRE_EQ(zs.size(), h.size() / 2);

    auto result = aberth_autocorr_mt(h, zs, options);
    auto niter = result.first;
    auto converged = result.second;

    MESSAGE("aberth_autocorr MT degree-100: niter=", niter, " converged=", converged);
    CHECK_LE(niter, options.max_iters);
    CHECK_GT(niter, 0);
}

// =====================================================================
// Bairstow stress tests
// =====================================================================

TEST_CASE("stress test pbairstow_even with high-degree polynomial") {
    const auto h = random_poly(100);

    Options options;
    options.tolerance = 1e-9;
    options.max_iters = 4000;

    auto vrs = initial_guess(h);
    REQUIRE_EQ(vrs.size(), h.size() / 2);

    auto result = pbairstow_even(h, vrs, options);
    auto niter = result.first;
    auto converged = result.second;

    MESSAGE("pbairstow_even degree-100: niter=", niter, " converged=", converged);
    CHECK_LE(niter, options.max_iters);
    CHECK_GT(niter, 0);
}

TEST_CASE("stress test pbairstow_autocorr with high-degree palindromic") {
    const auto h = random_palindromic_poly(100);

    Options options;
    options.tolerance = 1e-9;
    options.max_iters = 4000;

    auto vrs = initial_autocorr(h);
    REQUIRE_EQ(vrs.size(), h.size() / 4);

    auto result = pbairstow_autocorr(h, vrs, options);
    auto niter = result.first;
    auto converged = result.second;

    MESSAGE("pbairstow_autocorr degree-100: niter=", niter, " converged=", converged);
    CHECK_LE(niter, options.max_iters);
    CHECK_GT(niter, 0);
}

// =====================================================================
// Many small polynomial stress tests
// =====================================================================

TEST_CASE("stress test 100 random aberth polynomials") {
    Options options;
    options.tolerance = 1e-9;
    options.max_iters = 1000;

    int converged_count = 0;
    for (int trial = 0; trial < 100; ++trial) {
        auto h = random_poly(6);
        auto zs = initial_aberth(h);
        auto result = aberth(h, zs, options);
        if (result.second) {
            ++converged_count;
        }
    }
    // At least 80% should converge
    CHECK_GE(converged_count, 80);
}

TEST_CASE("stress test 100 random bairstow polynomials") {
    Options options;
    options.tolerance = 1e-9;
    options.max_iters = 1000;

    int converged_count = 0;
    for (int trial = 0; trial < 100; ++trial) {
        auto h = random_poly(6);
        auto vrs = initial_guess(h);
        auto result = pbairstow_even(h, vrs, options);
        if (result.second) {
            ++converged_count;
        }
    }
    // At least 80% should converge
    CHECK_GE(converged_count, 80);
}
