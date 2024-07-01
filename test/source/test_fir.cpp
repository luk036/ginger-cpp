// import numpy as np
// -*- coding: utf-8 -*-
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST_CASE

#include <ginger/autocorr.hpp>     // for extract_autocorr, initial_autocorr
#include <ginger/config.hpp>       // for Options
#include <ginger/rootfinding.hpp>  // for horner, Options
#include <utility>                 // for pair
#include <vector>                  // for vector

// #include "fmt/format.h"        // for print
// #include "ginger/vector2.hpp"  // for vector2

TEST_CASE("test FIR") {
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
    auto coeffs1 = r;
    // auto N = coeffs1.size() - 1;
    // auto vAh = horner(coeffs1, N, vrs[1]);
    // fmt::print("{}, {}\n", vAh.x(), vAh.y());
    // fmt::print(coeffs1);
    // auto vA1h = horner(coeffs1, N - 2, vrs[1]);
    // fmt::print("{}, {}\n", vA1h.x(), vA1h.y());

    auto options = Options();
    options.tolerance = 1e-5;
    auto result = pbairstow_autocorr(r, vrs, options);
    // auto niter = result.first;
    auto found = result.second;
    // fmt::print("{}, {}\n", niter, found);
    REQUIRE(found);

    // for (auto &vr : vrs) {
    //     // extract_autocorr(vr);
    //     fmt::print("{}, {}\n", vr.x(), vr.y());
    // }
}
