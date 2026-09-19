#include <doctest/doctest.h>  // for CHECK, REQUIRE, TEST_CASE

#include <algorithm>               // for sort, min
#include <cmath>                   // for abs
#include <complex>                 // for complex
#include <cstddef>                 // for size_t
#include <ginger/aberth.hpp>       // for aberth, initial_aberth, leja_order, poly_from_roots
#include <ginger/config.hpp>       // for ginger::Options, residual_scale
#include <ginger/rootfinding.hpp>  // for initial_guess, pbairstow_even
#include <limits>                  // for numeric_limits
#include <random>                  // for mt19937, uniform_real_distribution
#include <vector>                  // for vector

namespace {

    using Complex = std::complex<double>;

    auto real_poly(std::size_t degree, unsigned seed = 42) -> std::vector<double> {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> re(-3.0, 3.0);
        std::uniform_real_distribution<double> im(0.2, 3.0);
        std::vector<Complex> roots;
        roots.reserve(degree);
        while (roots.size() < degree) {
            const auto a = re(rng);
            const auto b = im(rng);
            roots.emplace_back(a, b);
            if (roots.size() < degree) {
                roots.emplace_back(a, -b);
            }
        }
        return poly_from_roots(roots);
    }

    auto random_points(std::size_t n, unsigned seed = 7) -> std::vector<Complex> {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> d(0.0, 1.0);
        std::vector<Complex> pts;
        pts.reserve(n);
        for (auto i = std::size_t{0}; i < n; ++i) {
            pts.emplace_back(d(rng), d(rng));
        }
        return pts;
    }

    /// Reference O(n^3) greedy Leja ordering used to validate leja_order.
    auto naive_leja(const std::vector<Complex>& points) -> std::vector<Complex> {
        if (points.empty()) {
            return {};
        }
        auto sorted = points;
        std::sort(sorted.begin(), sorted.end(),
                  [](const Complex& a, const Complex& b) { return std::abs(a) < std::abs(b); });
        std::vector<Complex> result;
        result.push_back(sorted.front());
        sorted.erase(sorted.begin());
        while (!sorted.empty()) {
            auto best_idx = std::size_t{0};
            auto best_dist = -1.0;
            for (auto i = std::size_t{0}; i < sorted.size(); ++i) {
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
            sorted.erase(sorted.begin() + static_cast<std::ptrdiff_t>(best_idx));
        }
        return result;
    }

    auto pre_scaled(const std::vector<double>& coeffs) -> std::vector<double> {
        const auto scale = ginger::residual_scale(coeffs);
        auto out = coeffs;
        for (auto& c : out) {
            c /= scale;
        }
        return out;
    }

}  // namespace

TEST_CASE("residual_scale is floored at one") {
    CHECK(ginger::residual_scale({0.1, 0.2}) == doctest::Approx(1.0));
    CHECK(ginger::residual_scale({10.0, 150.0, 10.0}) == doctest::Approx(150.0));
    CHECK(ginger::residual_scale({-3.0, 2.0}) == doctest::Approx(3.0));
}

TEST_CASE("leja_order matches the naive greedy reference") {
    for (const auto n : {2U, 3U, 10U, 37U, 100U}) {
        const auto pts = random_points(n);
        const auto fast = leja_order(pts);
        const auto slow = naive_leja(pts);
        REQUIRE(fast.size() == slow.size());
        for (auto i = std::size_t{0}; i < fast.size(); ++i) {
            CHECK(fast[i] == slow[i]);
        }
    }
}

TEST_CASE("leja_order returns a permutation of its input") {
    const auto pts = random_points(64);
    auto result = leja_order(pts);
    REQUIRE(result.size() == pts.size());
    const auto by_value = [](const Complex& a, const Complex& b) {
        return a.real() < b.real() || (a.real() == b.real() && a.imag() < b.imag());
    };
    auto sorted_result = result;
    auto sorted_input = pts;
    std::sort(sorted_result.begin(), sorted_result.end(), by_value);
    std::sort(sorted_input.begin(), sorted_input.end(), by_value);
    for (auto i = std::size_t{0}; i < sorted_input.size(); ++i) {
        CHECK(sorted_result[i] == sorted_input[i]);
    }
}

TEST_CASE("leja_order handles degenerate sizes") {
    CHECK(leja_order({}).empty());
    const auto one = leja_order({Complex{1.0, 2.0}});
    REQUIRE(one.size() == 1);
    CHECK(one[0] == Complex{1.0, 2.0});
}

TEST_CASE("aberth converges on a badly scaled polynomial at default tolerance") {
    const auto coeffs = real_poly(16);
    REQUIRE(ginger::residual_scale(coeffs) > 1e4);
    auto zs = initial_aberth(coeffs);
    const auto [niter, ok] = aberth(coeffs, zs, ginger::Options());
    CHECK(ok);
    CHECK(niter < 20);
}

TEST_CASE("coefficient scaling leaves the aberth roots unchanged") {
    const auto coeffs = real_poly(16);
    const auto zs0 = initial_aberth(coeffs);
    auto zs_a = zs0;
    auto zs_b = zs0;
    const auto [niter_a, ok_a] = aberth(coeffs, zs_a, ginger::Options());
    const auto [niter_b, ok_b] = aberth(pre_scaled(coeffs), zs_b, ginger::Options());
    CHECK(ok_a);
    CHECK(ok_b);
    CHECK(niter_a == niter_b);
    for (auto i = std::size_t{0}; i < zs_a.size(); ++i) {
        CHECK(std::abs(zs_a[i] - zs_b[i]) < 1e-9);
    }
}

TEST_CASE("coefficient scaling leaves the bairstow factors unchanged") {
    const auto coeffs = real_poly(8);
    const auto vrs0 = initial_guess(coeffs);
    auto vrs_a = vrs0;
    auto vrs_b = vrs0;
    pbairstow_even(coeffs, vrs_a, ginger::Options());
    pbairstow_even(pre_scaled(coeffs), vrs_b, ginger::Options());
    REQUIRE(vrs_a.size() == vrs_b.size());
    for (auto i = std::size_t{0}; i < vrs_a.size(); ++i) {
        CHECK(std::abs(vrs_a[i].x() - vrs_b[i].x()) < 1e-9);
        CHECK(std::abs(vrs_a[i].y() - vrs_b[i].y()) < 1e-9);
    }
}

TEST_CASE("bairstow converges on the degree-8 palindromic polynomial") {
    const auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_guess(h);
    const auto [niter, ok] = pbairstow_even(h, vrs, ginger::Options());
    CHECK(ok);
    CHECK(niter <= 12);
}
