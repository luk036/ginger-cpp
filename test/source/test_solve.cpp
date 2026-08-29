// Tests for the solve() facade entry points (auto policy selection).
#include <doctest/doctest.h>

#include <ginger/config.hpp>  // for ginger::Options, should_parallelize
#include <ginger/solve.hpp>   // for solve_aberth, solve_pbairstow_even, solve_mode

using namespace ginger;

TEST_CASE("test should_parallelize threshold") {
    CHECK_FALSE(should_parallelize(0));
    CHECK_FALSE(should_parallelize(4));
    CHECK(should_parallelize(5));
}

TEST_CASE("test solve facade pbairstow_even") {
    const auto h
        = std::vector<double>{1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 3.0, 0.0, 2.0, 0.0, 1.0};
    auto options = ginger::Options();
    options.tolerance = 1e-12;

    // 6 factors -> automatic dispatches to the multi-threaded variant.
    auto vrs = initial_guess(h);
    auto result = solve_pbairstow_even(h, vrs, options);
    REQUIRE(result.second);
    auto monic = poly_from_quadratic_factors(vrs);
    REQUIRE(monic.size() == h.size());
    auto scale = h[0];
    for (auto i = 0U; i < h.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(h[i]).epsilon(1e-8));
    }

    // Explicit mode selection must also converge.
    auto vrs2 = initial_guess(h);
    auto seq = solve_pbairstow_even(h, vrs2, options, solve_mode::sequential);
    CHECK(seq.second);
    auto vrs3 = initial_guess(h);
    auto atom = solve_pbairstow_even(h, vrs3, options, solve_mode::atomic);
    CHECK(atom.second);
}

TEST_CASE("test solve facade aberth") {
    const auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto options = ginger::Options();
    options.tolerance = 1e-12;

    // 4 roots -> automatic dispatches to the single-threaded variant.
    auto zs = initial_aberth(h);
    auto result = solve_aberth(h, zs, options);
    REQUIRE(result.second);
    auto monic = poly_from_roots(zs);
    REQUIRE(monic.size() == h.size());
    auto scale = h[0];
    for (auto i = 0U; i < h.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(h[i]).epsilon(1e-8));
    }

    auto zs2 = initial_aberth(h);
    auto mt = solve_aberth(h, zs2, options, solve_mode::multi_threaded);
    CHECK(mt.second);
    auto zs3 = initial_aberth(h);
    auto atom = solve_aberth(h, zs3, options, solve_mode::atomic);
    CHECK(atom.second);
}
