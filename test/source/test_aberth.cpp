// import numpy as np
// -*- coding: utf-8 -*-
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST_CASE

#include <ginger/aberth.hpp>       // for aberth, initial_aberth
#include <ginger/rootfinding.hpp>  // for Options
#include <utility>                 // for pair
#include <vector>                  // for vector

#include "fmt/format.h"  // for print

TEST_CASE("test aberth 1") {
    auto h = std::vector<double>{5., 2., 9., 6., 2.};
    auto zs = initial_aberth(h);
    auto result = aberth(h, zs, Options());
    auto niter = result.first;
    // auto found = result.second;
    // fmt::print("{}, {}\n", niter, found);

    CHECK(niter <= 11);
    // fmt::print([find_rootq(-r[0], -r[1]) for r : zs]);
}

TEST_CASE("test aberth 2") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto zs = initial_aberth(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth(h, zs, options);
    auto niter = result.first;
    // auto found = result.second;
    // fmt::print("{}, {}\n", niter, found);

    CHECK(niter <= 13);
    // fmt::print([find_rootq(-r[0], -r[1]) for r : zs]);
}

TEST_CASE("test aberth FIR") {
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
    auto result = aberth(r, zs, options);
    auto niter = result.first;
    // auto found = result.second;
    // fmt::print("{}, {}\n", niter, found);

    CHECK(niter <= 12);
    // fmt::print([find_rootq(-r[0], -r[1]) for r : zs]);
}

TEST_CASE("test horners method") {
    auto r = std::vector<double>{10, 34, 75, 94, 150, 94, 75, 34, 10};
    auto zs = initial_aberth(r);
    auto options = Options();
    options.tolerance = 1e-8;
    auto result = aberth(r, zs, options);
    auto niter = result.first;
    // auto found = result.second;
    // fmt::print("{}, {}\n", niter, found);
    CHECK(niter <= 12);
}
