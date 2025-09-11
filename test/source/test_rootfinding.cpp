// import numpy as np
// -*- coding: utf-8 -*-
#define DOCTEST_CONFIG_NO_EXCEPTIONS_BUT_WITH_ALL_ASSERTS
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST_CASE

#include <ginger/config.hpp>       // for Options
#include <ginger/rootfinding.hpp>  // for horner, initial_guess, pbairstow...
#include <utility>                 // for pair
#include <vector>                  // for vector

// #include "fmt/format.h"        // for print
#include "ginger/vector2.hpp"  // for Vector2

using namespace ginger;

TEST_CASE("test delta1()") {
    const auto vri = Vector2<double>(-2.0, 0.0);
    const auto vrj = Vector2<double>(4.0, 5.0);
    const auto vrk = Vector2<double>(3.0, 7.0);
    const auto vpj = vri - vrj;
    const auto vpk = vri - vrk;

    auto vA = Vector2<double>(3.0, 3.0);
    vA = delta(vA, vri, vpj);
    auto dr1 = delta(vA, vri, vpk);

    vA = Vector2<double>(3.0, 3.0);
    vA = delta(vA, vri, vpk);
    auto dr2 = delta(vA, vri, vpj);
    CHECK_EQ(dr1.dot(dr1), doctest::Approx(dr2.dot(dr2)));
}

TEST_CASE("test suppress 1") {
    const auto vri = Vector2<double>(-2.0, 0.0);
    const auto vrj = Vector2<double>(4.0, 5.0);
    const auto vrk = Vector2<double>(3.0, 7.0);

    auto vA = Vector2<double>(3.0, 3.0);
    auto vA1 = Vector2<double>(1.0, 2.0);

    suppress(vA, vA1, vri, vrj);
    suppress(vA, vA1, vri, vrk);
    auto dr1 = vA;

    vA = Vector2<double>(3.0, 3.0);
    vA1 = Vector2<double>(1.0, 2.0);
    suppress(vA, vA1, vri, vrk);
    suppress(vA, vA1, vri, vrj);
    auto dr2 = vA;
    CHECK_EQ(dr1.dot(dr1), doctest::Approx(dr2.dot(dr2)));
}

TEST_CASE("test suppress 2") {
    const auto vri = Vector2<double>(-2.0, 0.0);
    const auto vrj = Vector2<double>(4.0, 5.0);
    const auto vrk = Vector2<double>(3.0, 7.0);
    const auto vrl = Vector2<double>(3.0, 7.0);

    auto vA = Vector2<double>(3.0, 3.0);
    auto vA1 = Vector2<double>(1.0, 2.0);

    suppress(vA, vA1, vri, vrj);
    suppress(vA, vA1, vri, vrk);
    suppress(vA, vA1, vri, vrl);
    auto dr1 = vA;

    vA = Vector2<double>(3.0, 3.0);
    vA1 = Vector2<double>(1.0, 2.0);
    suppress(vA, vA1, vri, vrl);
    suppress(vA, vA1, vri, vrk);
    suppress(vA, vA1, vri, vrj);
    auto dr2 = vA;
    CHECK_EQ(dr1.dot(dr1), doctest::Approx(dr2.dot(dr2)));
}

TEST_CASE("test suppress 3") {
    const auto vri = Vector2<double>(-2.0, 0.0);
    const auto vrj = Vector2<double>(4.0, 5.0);
    const auto vrk = Vector2<double>(3.0, 7.0);

    auto vA = Vector2<double>(3.0, 3.0);
    auto vA1 = Vector2<double>(1.0, 2.0);

    suppress2(vA, vA1, vri, vrj);
    suppress2(vA, vA1, vri, vrk);
    auto dr1 = vA;

    vA = Vector2<double>(3.0, 3.0);
    vA1 = Vector2<double>(1.0, 2.0);
    suppress2(vA, vA1, vri, vrk);
    suppress2(vA, vA1, vri, vrj);
    auto dr2 = vA;
    CHECK_EQ(dr1.dot(dr1), doctest::Approx(dr2.dot(dr2)));
}

TEST_CASE("test suppress 4") {
    const auto vri = Vector2<double>(-2.0, 0.0);
    const auto vrj = Vector2<double>(4.0, 5.0);
    const auto vrk = Vector2<double>(3.0, 7.0);
    const auto vrl = Vector2<double>(3.0, 7.0);

    auto vA = Vector2<double>(3.0, 3.0);
    auto vA1 = Vector2<double>(1.0, 2.0);

    suppress2(vA, vA1, vri, vrj);
    suppress2(vA, vA1, vri, vrk);
    suppress2(vA, vA1, vri, vrl);
    auto dr1 = vA;

    vA = Vector2<double>(3.0, 3.0);
    vA1 = Vector2<double>(1.0, 2.0);
    suppress2(vA, vA1, vri, vrl);
    suppress2(vA, vA1, vri, vrk);
    suppress2(vA, vA1, vri, vrj);
    auto dr2 = vA;
    CHECK_EQ(dr1.dot(dr1), doctest::Approx(dr2.dot(dr2)));
}

TEST_CASE("test suppress 5") {
    const auto vri = Vector2<double>(-2.0, 0.0);
    const auto vrj = Vector2<double>(4.0, 5.0);
    // const auto vrk = Vector2<double>(3.0, 7.0);

    auto vA = Vector2<double>(3.0, 3.0);
    auto vA1 = Vector2<double>(1.0, 2.0);

    suppress(vA, vA1, vri, vrj);
    // suppress(vA, vA1, vri, vrk);
    auto dr1 = delta(vA, vri, vA1);

    vA = Vector2<double>(3.0, 3.0);
    vA1 = Vector2<double>(1.0, 2.0);
    suppress2(vA, vA1, vri, vrj);
    // suppress2(vA, vA1, vri, vrk);
    auto dr2 = delta(vA, vri, vA1);
    CHECK_EQ(dr1.dot(dr1), doctest::Approx(dr2.dot(dr2)));
}

TEST_CASE("test horner_eval") {
    auto h = std::vector<double>{1.0, 2.0, 3.0, 4.0};
    auto vA = horner_eval(h, 3, 1.0);
    CHECK(vA == 10.0);
    auto vA1 = horner_eval(h, 2, 1.0);
    CHECK(vA1 == 6.0);
}

TEST_CASE("test horner") {
    auto h = std::vector<double>{1.0, 2.0, 3.0, 4.0};
    auto v = ginger::Vector2<double>(1.0, 2.0);
    auto vA = horner(h, 3, v);
    CHECK(vA.x() == 8.0);
    CHECK(vA.y() == 10.0);
}

TEST_CASE("test initial_guess") {
    auto h = std::vector<double>{5., 2., 9., 6., 2.};
    auto vrs = initial_guess(h);
    CHECK(vrs.size() == 2);
}

TEST_CASE("test root-finding 1") {
    auto h = std::vector<double>{5., 2., 9., 6., 2.};
    auto vrs = initial_guess(h);
    // fmt::print(vrs);
    // fmt::print("vrs[1]: {}, {}\n", vrs[1].x(), vrs[1].y());
    auto coeffs1 = h;
    // auto degree = coeffs1.size() - 1;
    // auto vAh = horner(coeffs1, degree, vrs[1]);
    // fmt::print("{}, {}\n", vAh.x(), vAh.y());
    // fmt::print(coeffs1);
    // auto vA1h = horner(coeffs1, degree - 2, vrs[1]);
    // fmt::print("{}, {}\n", vA1h.x(), vA1h.y());

    auto result = pbairstow_even(h, vrs, Options());
    auto niter = result.first;
    auto found = result.second;
    // fmt::print("{}, {}\n", niter, found);

    REQUIRE(found);
    CHECK(niter <= 11);

    // fmt::print([find_rootq(-r[0], -r[1]) for r : vrs]);
}

TEST_CASE("test root-finding 2") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_guess(h);
    // fmt::print(vrs);
    // fmt::print("vrs[1]: {}, {}\n", vrs[1].x(), vrs[1].y());
    auto coeffs1 = h;
    // auto degree = coeffs1.size() - 1;
    // auto vAh = horner(coeffs1, degree, vrs[1]);
    // fmt::print("{}, {}\n", vAh.x(), vAh.y());
    // fmt::print(coeffs1);
    // auto vA1h = horner(coeffs1, degree - 2, vrs[1]);
    // fmt::print("{}, {}\n", vA1h.x(), vA1h.y());

    auto options = Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_even(h, vrs, options);
    auto niter = result.first;
    auto found = result.second;
    // fmt::print("{}, {}\n", niter, found);

    REQUIRE(found);
    CHECK(niter <= 13);
    // fmt::print([find_rootq(-r[0], -r[1]) for r : vrs]);
}

TEST_CASE("test root-finding FIR") {
    // auto vA = vec2{0.1, 1.2};
    // auto vA1 = vec2{2.3, 3.4};
    // auto vr = vec2{4.5, 5.6};
    // auto vrj = vec2{6.7, 7.8};
    // auto vA1 = suppress(vA, vA1, vr, vrj);
    // fmt::print(check_newton(vA, vA1, vr));
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_guess(h);
    // fmt::print(vrs);
    // fmt::print("vrs[1]: {}, {}\n", vrs[1].x(), vrs[1].y());
    auto coeffs1 = h;
    // auto degree = coeffs1.size() - 1;
    // auto vAh = horner(coeffs1, degree, vrs[1]);
    // fmt::print("{}, {}\n", vAh.x(), vAh.y());
    // fmt::print(coeffs1);
    // auto vA1h = horner(coeffs1, degree - 2, vrs[1]);
    // fmt::print("{}, {}\n", vA1h.x(), vA1h.y());

    auto options = Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_even(h, vrs, options);
    auto niter = result.first;
    // auto found = result.second;
    // fmt::print("{}, {}\n", niter, found);

    CHECK(niter <= 14);
    // fmt::print([find_rootq(-r[0], -r[1]) for r : vrs]);
}

TEST_CASE("Polynomial Root Finding") {
    SUBCASE("Initial Guess") {
        std::vector<double> coeffs = {10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
        auto vrs = initial_guess(coeffs);

        CHECK(vrs.size() == 4);  // For degree 8 polynomial
        for (const auto& vr : vrs) {
            CHECK(vr.x() != 0.0);
            CHECK(vr.y() != 0.0);
        }
    }

    SUBCASE("Horner Evaluation") {
        std::vector<double> coeffs = {1.0, 2.0, 3.0};  // x^2 + 2x + 3
        double z = 2.0;
        double result = horner_eval(coeffs, 2, z);
        CHECK(result == doctest::Approx(z * z + 2.0 * z + 3.0));
    }

    SUBCASE("Matrix Operations") {
        Vec2 vr(1.0, 2.0);
        Vec2 vp(3.0, 4.0);

        SUBCASE("Make Adjoint") {
            Mat2 adj = makeadjoint(vr, vp);
            CHECK(adj.x().x() == doctest::Approx(4.0));
            CHECK(adj.x().y() == doctest::Approx(-3.0));
            CHECK(adj.y().x() == doctest::Approx(-6.0));
            CHECK(adj.y().y() == doctest::Approx(7.0));
        }

        SUBCASE("Delta Calculation") {
            Vec2 vA(1.0, 2.0);
            Vec2 dt = delta(vA, vr, vp);
            CHECK(dt.x() != 0.0);
            CHECK(dt.y() != 0.0);
        }
    }

    SUBCASE("Suppression Functions") {
        Vec2 vA(3.0, 3.0);
        Vec2 vA1(1.0, 2.0);
        Vec2 vri(-2.0, 0.0);
        Vec2 vrj(4.0, 5.0);

        SUBCASE("Original Suppress") {
            Vec2 vA_copy = vA;
            Vec2 vA1_copy = vA1;
            suppress(vA_copy, vA1_copy, vri, vrj);
            CHECK(vA_copy.x() != vA.x());
            CHECK(vA_copy.y() != vA.y());
            CHECK(vA1_copy.x() != vA1.x());
            CHECK(vA1_copy.y() != vA1.y());
        }

        SUBCASE("Suppress2") {
            Vec2 vA_copy = vA;
            Vec2 vA1_copy = vA1;
            suppress2(vA_copy, vA1_copy, vri, vrj);
            CHECK(vA_copy.x() != vA.x());
            CHECK(vA_copy.y() != vA.y());
            CHECK(vA1_copy.x() != vA1.x());
            CHECK(vA1_copy.y() != vA1.y());
        }
    }

    SUBCASE("Bairstow's Method") {
        std::vector<double> coeffs = {10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
        auto vrs = initial_guess(coeffs);
        Options options;
        options.max_iters = 100;
        options.tolerance = 1e-12;

        SUBCASE("Convergence") {
            auto [niter, found] = pbairstow_even(coeffs, vrs, options);
            CHECK(niter > 0);
            CHECK(found == true);

            // Verify roots by evaluating polynomial
            for (const auto& vr : vrs) {
                std::vector<double> coeffs_copy = coeffs;
                auto val = horner(coeffs_copy, coeffs.size() - 1, vr);
                CHECK(val.x() == doctest::Approx(0.0).epsilon(1e-6));
                CHECK(val.y() == doctest::Approx(0.0).epsilon(1e-6));
            }
        }

        SUBCASE("Max Iterations") {
            options.max_iters = 1;
            auto [niter, found] = pbairstow_even(coeffs, vrs, options);
            CHECK(niter == 1);
            CHECK(found == false);
        }
    }
}

// TEST_CASE("Edge Cases") {
//     SUBCASE("Zero Polynomial") {
//         std::vector<double> coeffs = {0.0, 0.0, 0.0};
//         auto vrs = initial_guess(coeffs);
//         CHECK(vrs.empty()); // Should handle zero polynomial
//
//         Options options;
//         auto [niter, found] = pbairstow_even(coeffs, vrs, options);
//         CHECK(found == true); // Technically correct for zero polynomial
//     }
//
//     SUBCASE("Constant Polynomial") {
//         std::vector<double> coeffs = {5.0};
//         auto vrs = initial_guess(coeffs);
//         CHECK(vrs.empty()); // No roots to find
//
//         Options options;
//         auto [niter, found] = pbairstow_even(coeffs, vrs, options);
//         CHECK(found == true); // No iterations needed
//     }
// }
