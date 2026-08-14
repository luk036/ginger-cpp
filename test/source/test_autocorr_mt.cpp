// import numpy as np
// -*- coding: utf-8 -*-
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST_CASE

#include <ginger/autocorr_mt.hpp>  // for extract_autocorr, initial_autocorr, pbairstow_...
#include <ginger/config.hpp>       // for Options
#include <ginger/vector2.hpp>      // for vector2
#include <utility>                 // for pair
#include <vector>                  // for vector

using namespace ginger;

TEST_CASE("test auto-corr mt 1") {
    auto r = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_autocorr(r);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_autocorr_mt(r, vrs, options);
    auto niter = result.first;
    auto found = result.second;
    for (auto& vr : vrs) {
        extract_autocorr(vr);
    }
    REQUIRE(found);
    CHECK_LE(niter, 21);
}

TEST_CASE("test autocorr mt FIR") {
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
    auto result = pbairstow_autocorr_mt(r, vrs, options);
    auto found = result.second;
    REQUIRE(found);
}

TEST_CASE("test poly_from_autocorr_factors mt reconstruction") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_autocorr(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_autocorr_mt(h, vrs, options);
    REQUIRE(result.second);
    auto monic = poly_from_autocorr_factors(vrs);
    REQUIRE(monic.size() == h.size());
    auto scale = h[0];
    for (auto i = 0U; i < h.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(h[i]).epsilon(1e-8));
    }
}

TEST_CASE("test poly_from_autocorr_factors mt fir") {
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
    auto result = pbairstow_autocorr_mt(r, vrs, options);
    REQUIRE(result.second);
    auto monic = poly_from_autocorr_factors(vrs);
    REQUIRE(monic.size() == r.size());
    auto scale = r[0];
    for (auto i = 0U; i < r.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(r[i]).epsilon(1e-3));
    }
}
