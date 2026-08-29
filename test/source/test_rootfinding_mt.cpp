// import numpy as np
// -*- coding: utf-8 -*-
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST_CASE

#include <ginger/config.hpp>          // for ginger::Options
#include <ginger/rootfinding_mt.hpp>  // for initial_guess, pbairstow_even_mt, poly_from_...
#include <utility>                    // for pair
#include <vector>                     // for vector

using namespace ginger;

TEST_CASE("test root-finding mt 1") {
    auto h = std::vector<double>{5., 2., 9., 6., 2.};
    auto vrs = initial_guess(h);
    auto result = pbairstow_even_mt(h, vrs, ginger::Options());
    auto niter = result.first;
    CHECK_LE(niter, 11);
}

TEST_CASE("test root-finding mt 2") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_guess(h);
    auto options = ginger::Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_even_mt(h, vrs, options);
    auto niter = result.first;
    CHECK_LE(niter, 13);
}

TEST_CASE("test root-finding mt FIR") {
    auto r = std::vector<double>{
        -0.00196191, -0.00094597, -0.00023823, 0.00134667,  0.00380494,  0.00681596,  0.0097864,
        0.01186197,  0.0121238,   0.00985211,  0.00474894,  -0.00281751, -0.01173923, -0.0201885,
        -0.02590168, -0.02658216, -0.02035729, -0.00628271, 0.01534627,  0.04279982,  0.0732094,
        0.10275561,  0.12753013,  0.14399228,  0.15265722,  0.14399228,  0.12753013,  0.10275561,
        0.0732094,   0.04279982,  0.01534627,  -0.00628271, -0.02035729, -0.02658216, -0.02590168,
        -0.0201885,  -0.01173923, -0.00281751, 0.00474894,  0.00985211,  0.0121238,   0.01186197,
        0.0097864,   0.00681596,  0.00380494,  0.00134667,  -0.00023823, -0.00094597, -0.00196191,
    };
    auto vrs = initial_guess(r);
    auto options = ginger::Options();
    options.tolerance = 1e-6;
    auto result = pbairstow_even_mt(r, vrs, options);
    auto niter = result.first;
    CHECK_LE(niter, 14);
}

TEST_CASE("test root-finding mt degree 12") {
    // Palindromic degree-12 polynomial: 6 quadratic factors, exercises the
    // true multi-threaded path (num_roots > 4).
    auto h = std::vector<double>{1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 3.0, 0.0, 2.0, 0.0, 1.0};
    auto vrs = initial_guess(h);
    auto options = ginger::Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_even_mt(h, vrs, options);
    auto niter = result.first;
    CHECK_LE(niter, 14);
}

TEST_CASE("test poly_from_quadratic_factors mt reconstruction") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_guess(h);
    auto options = ginger::Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_even_mt(h, vrs, options);
    REQUIRE(result.second);
    auto monic = poly_from_quadratic_factors(vrs);
    REQUIRE(monic.size() == h.size());
    auto scale = h[0];
    for (auto i = 0U; i < h.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(h[i]).epsilon(1e-8));
    }
}
