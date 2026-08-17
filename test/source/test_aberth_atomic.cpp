// import numpy as np
// -*- coding: utf-8 -*-
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST_CASE

#include <ginger/aberth_atomic.hpp>  // for aberth_atomic, initial_aberth
#include <ginger/config.hpp>         // for Options
#include <utility>                   // for pair
#include <vector>                    // for vector

TEST_CASE("test aberth_atomic 1") {
    auto h = std::vector<double>{5., 2., 9., 6., 2.};
    auto zs = initial_aberth(h);
    auto result = aberth_atomic(h, zs, Options());
    auto niter = result.first;
    CHECK(result.second);
    CHECK_LE(niter, 12);
}

TEST_CASE("test aberth_atomic 2") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto zs = initial_aberth(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth_atomic(h, zs, options);
    auto niter = result.first;
    CHECK(result.second);
    CHECK_LE(niter, 30);
}

TEST_CASE("test aberth_atomic FIR") {
    auto r = std::vector<double>{
        -0.00196191, -0.00094597, -0.00023823, 0.00134667,  0.00380494,  0.00681596,  0.0097864,
        0.01186197,  0.0121238,   0.00985211,  0.00474894,  -0.00281751, -0.01173923, -0.0201885,
        -0.02590168, -0.02658216, -0.02035729, -0.00628271, 0.01534627,  0.04279982,  0.0732094,
        0.10275561,  0.12753013,  0.14399228,  0.15265722,  0.14399228,  0.12753013,  0.10275561,
        0.0732094,   0.04279982,  0.01534627,  -0.00628271, -0.02035729, -0.02658216, -0.02590168,
        -0.0201885,  -0.01173923, -0.00281751, 0.00474894,  0.00985211,  0.0121238,   0.01186197,
        0.0097864,   0.00681596,  0.00380494,  0.00134667,  -0.00023823, -0.00094597, -0.00196191,
    };
    auto zs = initial_aberth(r);
    auto options = Options();
    options.tolerance = 1e-8;
    auto result = aberth_atomic(r, zs, options);
    auto niter = result.first;
    CHECK(result.second);
    CHECK_LE(niter, 40);
}

TEST_CASE("test aberth_autocorr_atomic 2") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto zs = initial_aberth_autocorr(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth_autocorr_atomic(h, zs, options);
    auto niter = result.first;
    CHECK(result.second);
    CHECK_LE(niter, 14);
}

TEST_CASE("test aberth_autocorr_atomic FIR") {
    auto r = std::vector<double>{
        -0.00196191, -0.00094597, -0.00023823, 0.00134667,  0.00380494,  0.00681596,  0.0097864,
        0.01186197,  0.0121238,   0.00985211,  0.00474894,  -0.00281751, -0.01173923, -0.0201885,
        -0.02590168, -0.02658216, -0.02035729, -0.00628271, 0.01534627,  0.04279982,  0.0732094,
        0.10275561,  0.12753013,  0.14399228,  0.15265722,  0.14399228,  0.12753013,  0.10275561,
        0.0732094,   0.04279982,  0.01534627,  -0.00628271, -0.02035729, -0.02658216, -0.02590168,
        -0.0201885,  -0.01173923, -0.00281751, 0.00474894,  0.00985211,  0.0121238,   0.01186197,
        0.0097864,   0.00681596,  0.00380494,  0.00134667,  -0.00023823, -0.00094597, -0.00196191,
    };
    auto zs = initial_aberth_autocorr(r);
    auto options = Options();
    options.tolerance = 1e-8;
    auto result = aberth_autocorr_atomic(r, zs, options);
    auto niter = result.first;
    CHECK(result.second);
    CHECK_LE(niter, 40);
}

TEST_CASE("test poly_from_autocorr_roots atomic reconstruction") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto zs = initial_aberth_autocorr(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth_autocorr_atomic(h, zs, options);
    REQUIRE(result.second);
    auto monic = poly_from_autocorr_roots(zs);
    REQUIRE(monic.size() == h.size());
    auto scale = h[0];
    for (auto i = 0U; i < h.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(h[i]).epsilon(1e-8));
    }
}
