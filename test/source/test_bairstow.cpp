// import numpy as np
// -*- coding: utf-8 -*-
#define DOCTEST_CONFIG_NO_EXCEPTIONS_BUT_WITH_ALL_ASSERTS
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST_CASE

#include <ginger/bairstow.hpp>     // for horner, initial_guess, pbairstow...
#include <ginger/config.hpp>       // for Options
#include <ginger/rootfinding.hpp>  // for horner, initial_guess, pbairstow...
// #include <utility>                 // for pair
#include <vector>  // for vector

// #include "fmt/format.h"        // for print
#include "ginger/vector2.hpp"  // for Vector2

using namespace ginger;

// TEST_CASE("test delta_ref") {
//     double x1 = 3;
//     double y1 = 3;
//     const auto vri = Vector2<double>(-2.0, 0.0);
//     const auto vrj = Vector2<double>(4.0, 5.0);
//     const auto vrk = Vector2<double>(3.0, 7.0);
//     const auto vpj = vri - vrj;
//     const auto vpk = vri - vrk;

//     auto vA = Vec2Ref(x1, y1);
//     vA = delta_ref(vA, vri, vpj);
//     auto dr1 = delta_ref(vA, vri, vpk);

//     double x2 = 3;
//     double y2 = 3;
//     auto vA2 = Vec2Ref(x2, y2);
//     vA2 = delta_ref(vA2, vri, vpk);
//     auto dr2 = delta_ref(vA, vri, vpj);
//     CHECK_EQ(dr1.dot(dr1), doctest::Approx(dr2.dot(dr2)));
// }

TEST_CASE("test horner_rule 1") {
    const auto v = ginger::Vector2<double>(1.0, 2.0);
    const auto coeffs = std::vector<double>{1.0, 2.0, 3.0, 4.0};
    const auto degree = coeffs.size() - 1;

    auto coeffs1 = coeffs;
    std::vector<Vec2Ref> vcoeffs;
    for (auto i = 0U; i < degree; ++i) {
        vcoeffs.emplace_back(Vec2Ref{coeffs1[i], coeffs1[i + 1]});
    }

    auto vA = horner_ref(coeffs1, vcoeffs, degree, v);
    CHECK_EQ(vA.x(), 8.0);
    CHECK_EQ(vA.y(), 10.0);

    auto coeffs2 = coeffs;
    auto vA2 = horner(coeffs2, degree, v);
    CHECK_EQ(vA2.x(), 8.0);
    CHECK_EQ(vA2.y(), 10.0);

    CHECK_EQ(coeffs1, coeffs2);
}

TEST_CASE("test horner_rule 2") {
    const auto vr = ginger::Vector2<double>(1.0, 2.0);
    const auto coeffs = std::vector<double>{5.0, 2.0, 9.0, 6.0, 2.0};
    const auto degree = coeffs.size() - 1;

    auto coeffs1 = coeffs;
    std::vector<Vec2Ref> vcoeffs;
    for (auto i = 0U; i < degree; ++i) {
        vcoeffs.emplace_back(Vec2Ref{coeffs1[i], coeffs1[i + 1]});
    }
    auto vA = horner_ref(coeffs1, vcoeffs, degree, vr);
    auto vA1 = horner_ref(coeffs1, vcoeffs, degree - 2, vr);
    auto vd1 = delta_ref(vA, vr, vA1);

    auto coeffs2 = coeffs;
    auto vA2 = horner(coeffs2, degree, vr);
    auto vA3 = horner(coeffs2, degree - 2, vr);
    auto vd2 = delta(vA2, vr, vA3);

    CHECK_EQ(vA.x(), vA2.x());
    CHECK_EQ(vA.y(), vA2.y());
    CHECK_EQ(vA1.x(), vA3.x());
    CHECK_EQ(vA1.y(), vA3.y());
    CHECK_EQ(coeffs1, coeffs2);
    CHECK_EQ(vd1.x(), vd2.x());
    CHECK_EQ(vd1.y(), vd2.y());
}

// TEST_CASE("test bairstow 1") {
//     auto h = std::vector<double>{5.0, 2.0, 9.0, 6.0, 2.0};
//     auto vr = Vector2{1.1, 2.2};
//     auto options = Options{2000, 1e-8};
//     auto result = bairstow(h, vr, options);
//     auto niter = result.first;
//     auto found = result.second;
//     fmt::print("-------------------- {}, {}\n", niter, found);
//     REQUIRE(found);
//     CHECK(niter <= 11);
//     // fmt::print([find_rootq(-r[0], -r[1]) for r : vrs]);
// }

// TEST_CASE("test bairstow 2") {
//     auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
//     auto vrs = initial_guess(h);
//     // fmt::print(vrs);
//     fmt::print("vrs[1]: {}, {}\n", vrs[1].x(), vrs[1].y());
//     auto coeffs1 = h;
//     auto degree = coeffs1.size() - 1;
//     auto vAh = horner(coeffs1, degree, vrs[1]);
//     fmt::print("{}, {}\n", vAh.x(), vAh.y());
//     // fmt::print(coeffs1);
//     auto vA1h = horner(coeffs1, degree - 2, vrs[1]);
//     fmt::print("{}, {}\n", vA1h.x(), vA1h.y());

//     auto options = Options();
//     options.tolerance = 1e-12;
//     auto result = pbairstow_even(h, vrs, options);
//     auto niter = result.first;
//     auto found = result.second;
//     fmt::print("{}, {}\n", niter, found);

//     REQUIRE(found);
//     CHECK(niter <= 13);
//     // fmt::print([find_rootq(-r[0], -r[1]) for r : vrs]);
// }

// TEST_CASE("test root-finding FIR") {
//     // auto vA = vec2{0.1, 1.2};
//     // auto vA1 = vec2{2.3, 3.4};
//     // auto vr = vec2{4.5, 5.6};
//     // auto vrj = vec2{6.7, 7.8};
//     // auto vA1 = suppress(vA, vA1, vr, vrj);
//     // fmt::print(check_newton(vA, vA1, vr));
//     auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
//     auto vrs = initial_guess(h);
//     // fmt::print(vrs);
//     fmt::print("vrs[1]: {}, {}\n", vrs[1].x(), vrs[1].y());
//     auto coeffs1 = h;
//     auto degree = coeffs1.size() - 1;
//     auto vAh = horner(coeffs1, degree, vrs[1]);
//     fmt::print("{}, {}\n", vAh.x(), vAh.y());
//     // fmt::print(coeffs1);
//     auto vA1h = horner(coeffs1, degree - 2, vrs[1]);
//     fmt::print("{}, {}\n", vA1h.x(), vA1h.y());

//     auto options = Options();
//     options.tolerance = 1e-12;
//     auto result = pbairstow_even(h, vrs, options);
//     auto niter = result.first;
//     auto found = result.second;
//     fmt::print("{}, {}\n", niter, found);

//     CHECK(niter <= 14);
//     // fmt::print([find_rootq(-r[0], -r[1]) for r : vrs]);
// }
