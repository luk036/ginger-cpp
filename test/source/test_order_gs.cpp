// Order-dependence tests for the single-threaded (Gauss-Seidel) Bairstow
// variant, plus the suppression-order drift test.
// Mirrors the Rust ginger-rs suite:
//   - pbairstow_even (Gauss-Seidel) is order-DEPENDENT in iteration count
//   - suppression order within a job only causes machine-epsilon drift
#include <doctest/doctest.h>

#include <algorithm>          // for all_of
#include <cmath>              // for abs
#include <ginger/config.hpp>  // for Options
#include <ginger/rootfinding.hpp>  // for delta_scalar, horner, initial_guess, pbairstow_even, suppress_old
#include <vector>                  // for vector

using namespace ginger;

TEST_CASE("test gs order dependent iterations") {
    const auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto options = Options();
    options.tolerance = 1e-12;

    const auto base = initial_guess(h);
    REQUIRE_EQ(base.size(), 4U);

    std::vector<std::vector<Vec2>> perms{
        base,
        {base[3], base[2], base[1], base[0]},
        {base[1], base[3], base[0], base[2]},
        {base[2], base[0], base[3], base[1]},
    };

    std::vector<unsigned int> niters;
    for (auto& perm : perms) {
        auto result = pbairstow_even(h, perm, options);
        REQUIRE(result.second);
        niters.push_back(result.first);
    }

    // Gauss-Seidel is NOT order-independent: iteration counts differ.
    const auto all_equal
        = std::all_of(niters.begin() + 1, niters.end(), [&](auto n) { return n == niters[0]; });
    CHECK_FALSE(all_equal);
}

TEST_CASE("test suppression order machine epsilon") {
    const auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    const auto vrs = initial_guess(h);

    std::vector<std::vector<size_t>> orders{
        {0, 1, 2, 3},
        {3, 2, 1, 0},
        {1, 3, 0, 2},
        {2, 0, 3, 1},
    };

    auto worst = 0.0;
    for (auto i = 0U; i < vrs.size(); ++i) {
        bool has_ref = false;
        Vec2 ref{0.0, 0.0};
        for (const auto& order : orders) {
            auto local_coeffs = h;  // horner corrupts the array
            const auto degree = local_coeffs.size() - 1;
            const auto& vri = vrs[i];
            auto vA = horner(local_coeffs, degree, vri);
            if (std::max(std::abs(vA.x()), std::abs(vA.y())) < 1e-15) {
                continue;
            }
            auto vA1 = horner(local_coeffs, degree - 2, vri);
            for (const auto j : order) {
                if (j != i) {
                    suppress_old(vA, vA1, vri, vrs[j]);
                }
            }
            const auto new_vri = vri - delta_scalar(vA, vri, vA1);
            if (!has_ref) {
                ref = new_vri;
                has_ref = true;
            } else {
                worst = std::max(worst, std::abs(new_vri.x() - ref.x()));
                worst = std::max(worst, std::abs(new_vri.y() - ref.y()));
            }
        }
    }
    CHECK_LT(worst, 1e-12);
}
