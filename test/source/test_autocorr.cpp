// import numpy as np
// -*- coding: utf-8 -*-
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST_CASE
#include <fmt/format.h>       // for print

#include <ginger/autocorr.hpp>     // for extract_autocorr, initial_autocorr
#include <ginger/config.hpp>       // for Options
#include <ginger/rootfinding.hpp>  // for horner, Options
#include <ginger/vector2.hpp>      // for vector2
#include <utility>                 // for pair
#include <vector>                  // for vector

TEST_CASE("test auto-corr 1") {
    // auto vA = vec2{0.1, 1.2};
    // auto vA1 = vec2{2.3, 3.4};
    // auto vr = vec2{4.5, 5.6};
    // auto vrj = vec2{6.7, 7.8};
    // auto vA1 = suppress(vA, vA1, vr, vrj);
    // fmt::print(check_newton(vA, vA1, vr));
    auto r = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_autocorr(r);
    // fmt::print(vrs);
    // fmt::print("vrs[1]: {}, {}\n", vrs[1].x(), vrs[1].y());
    // auto degree = coeffs1.size() - 1;
    // auto vAh = horner(coeffs1, degree, vrs[1]);
    // fmt::print("{}, {}\n", vAh.x(), vAh.y());
    // fmt::print(coeffs1);
    // auto vA1h = horner(coeffs1, degree - 2, vrs[1]);
    // fmt::print("{}, {}\n", vA1h.x(), vA1h.y());

    auto options = Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_autocorr(r, vrs, options);
    auto niter = result.first;
    auto found = result.second;
    // fmt::print("{}, {}\n", niter, found);
    for (auto& vr : vrs) {
        extract_autocorr(vr);
        fmt::print("{}, {}\n", vr.x(), vr.y());
    }
    REQUIRE(found);

    CHECK_LE(niter, 21);

    // fmt::print([find_rootq(-r[0], -r[1]) for r : vrs]);
}

TEST_CASE("test autocorr FIR") {
    // auto vA = vec2{0.1, 1.2};
    // auto vA1 = vec2{2.3, 3.4};
    // auto vr = vec2{4.5, 5.6};
    // auto vrj = vec2{6.7, 7.8};
    // auto vA1 = suppress(vA, vA1, vr, vrj);
    // fmt::print(check_newton(vA, vA1, vr));
    auto r = std::vector<double>{
        -0.00196191, -0.00094597, -0.00023823, 0.00134667,  0.00380494,  0.00681596,  0.0097864,
        0.01186197,  0.0121238,   0.00985211,  0.00474894,  -0.00281751, -0.01173923, -0.0201885,
        -0.02590168, -0.02658216, -0.02035729, -0.00628271, 0.01534627,  0.04279982,  0.0732094,
        0.10275561,  0.12753013,  0.14399228,  0.15265722,  0.14399228,  0.12753013,  0.10275561,
        0.0732094,   0.04279982,  0.01534627,  -0.00628271, -0.02035729, -0.02658216, -0.02590168,
        -0.0201885,  -0.01173923, -0.00281751, 0.00474894,  0.00985211,  0.0121238,   0.01186197,
        0.0097864,   0.00681596,  0.00380494,  0.00134667,  -0.00023823, -0.00094597, -0.00196191};
    auto vrs = initial_autocorr(r);
    // fmt::print(vrs);
    // fmt::print("vrs[1]: {}, {}\n", vrs[1].x(), vrs[1].y());
    // auto coeffs1 = r;
    // auto degree = coeffs1.size() - 1;
    // auto vAh = horner(coeffs1, degree, vrs[1]);
    // fmt::print("{}, {}\n", vAh.x(), vAh.y());
    // fmt::print(coeffs1);
    // auto vA1h = horner(coeffs1, degree - 2, vrs[1]);
    // fmt::print("{}, {}\n", vA1h.x(), vA1h.y());

    auto options = Options();
    options.tolerance = 1e-4;
    auto result = pbairstow_autocorr(r, vrs, options);
    // auto niter = result.first;
    auto found = result.second;
    // fmt::print("{}, {}\n", niter, found);
    REQUIRE(found);

    // for (auto &vr : vrs) {
    //     // extract_autocorr(vr);
    //     fmt::print("{}, {}\n", vr.x(), vr.y());
    // }

    // CHECK_LE(niter, 346);

    // fmt::print([find_rootq(-r[0], -r[1]) for r : vrs]);
}

TEST_CASE("test poly_from_autocorr_factors empty") {
    auto coeffs = poly_from_autocorr_factors({});
    REQUIRE(coeffs.size() == 1);
    CHECK_EQ(coeffs[0], doctest::Approx(1.0));
}

TEST_CASE("test poly_from_autocorr_factors reconstruction") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_autocorr(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_autocorr(h, vrs, options);
    REQUIRE(result.second);
    auto monic = poly_from_autocorr_factors(vrs);
    REQUIRE(monic.size() == h.size());
    auto scale = h[0];
    for (auto i = 0U; i < h.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(h[i]).epsilon(1e-8));
    }
}

TEST_CASE("test poly_from_autocorr_factors fir") {
    auto r = std::vector<double>{
        -0.00196191, -0.00094597, -0.00023823, 0.00134667,  0.00380494,  0.00681596,  0.0097864,
        0.01186197,  0.0121238,   0.00985211,  0.00474894,  -0.00281751, -0.01173923, -0.0201885,
        -0.02590168, -0.02658216, -0.02035729, -0.00628271, 0.01534627,  0.04279982,  0.0732094,
        0.10275561,  0.12753013,  0.14399228,  0.15265722,  0.14399228,  0.12753013,  0.10275561,
        0.0732094,   0.04279982,  0.01534627,  -0.00628271, -0.02035729, -0.02658216, -0.02590168,
        -0.0201885,  -0.01173923, -0.00281751, 0.00474894,  0.00985211,  0.0121238,   0.01186197,
        0.0097864,   0.00681596,  0.00380494,  0.00134667,  -0.00023823, -0.00094597, -0.00196191};
    auto vrs = initial_autocorr(r);
    auto options = Options();
    options.tolerance = 1e-4;
    auto result = pbairstow_autocorr(r, vrs, options);
    REQUIRE(result.second);
    auto monic = poly_from_autocorr_factors(vrs);
    REQUIRE(monic.size() == r.size());
    auto scale = r[0];
    for (auto i = 0U; i < r.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(r[i]).epsilon(1e-3));
    }
}
