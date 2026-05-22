// import numpy as np
// -*- coding: utf-8 -*-
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST_CASE

#include <ginger/aberth.hpp>  // for aberth, initial_aberth
#include <ginger/config.hpp>  // for Options
#include <utility>            // for pair
#include <vector>             // for vector

// #include "fmt/format.h"  // for print

TEST_CASE("test aberth 1") {
    auto h = std::vector<double>{5., 2., 9., 6., 2.};
    auto zs = initial_aberth(h);
    auto result = aberth(h, zs, Options());
    auto niter = result.first;
    CHECK_LE(niter, 11);
}

TEST_CASE("test aberth 2") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto zs = initial_aberth(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth(h, zs, options);
    auto niter = result.first;
    CHECK_LE(niter, 13);
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
    CHECK_LE(niter, 12);
}

TEST_CASE("test aberth_mt 1") {
    auto h = std::vector<double>{5., 2., 9., 6., 2.};
    auto zs = initial_aberth(h);
    auto result = aberth_mt(h, zs, Options());
    auto niter = result.first;
    CHECK_LE(niter, 12);
}

TEST_CASE("test aberth_mt 2") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto zs = initial_aberth(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth_mt(h, zs, options);
    auto niter = result.first;
    CHECK_LE(niter, 14);
}

TEST_CASE("test aberth_mt FIR") {
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
    auto result = aberth_mt(r, zs, options);
    auto niter = result.first;
    CHECK_LE(niter, 13);
}

// TEST_CASE("test aberth_autocorr 1") {
//     auto h = std::vector<double>{5., 2., 9., 6., 2.};
//     auto zs = initial_aberth_autocorr(h);
//     auto result = aberth_autocorr(h, zs, Options());
//     auto niter = result.first;
//     CHECK_LE(niter, 11);
// }

TEST_CASE("test aberth_autocorr 2") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto zs = initial_aberth_autocorr(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth(h, zs, options);
    auto niter = result.first;
    CHECK_LE(niter, 13);
}

TEST_CASE("test aberth_autocorr FIR") {
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
    auto result = aberth_autocorr(r, zs, options);
    auto niter = result.first;
    CHECK_LE(niter, 12);
}

// TEST_CASE("test aberth_autocorr_mt 1") {
//     auto h = std::vector<double>{5., 2., 9., 6., 2.};
//     auto zs = initial_aberth_autocorr(h);
//     auto result = aberth_autocorr_mt(h, zs, Options());
//     auto niter = result.first;
//     CHECK_LE(niter, 12);
// }

TEST_CASE("test aberth_autocorr_mt 2") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto zs = initial_aberth_autocorr(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth_mt(h, zs, options);
    auto niter = result.first;
    CHECK_LE(niter, 14);
}

TEST_CASE("test aberth_autocorr_mt FIR") {
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
    auto result = aberth_autocorr_mt(r, zs, options);
    auto niter = result.first;
    CHECK_LE(niter, 13);
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
    CHECK_LE(niter, 12);
}

TEST_CASE("test poly_from_roots known roots") {
    auto roots = std::vector<std::complex<double>>{{1.0, 0.0}, {-1.0, 0.0}};
    auto coeffs = poly_from_roots(roots);
    REQUIRE(coeffs.size() == 3);
    CHECK_EQ(coeffs[0], doctest::Approx(1.0));
    CHECK_EQ(coeffs[1], doctest::Approx(0.0));
    CHECK_EQ(coeffs[2], doctest::Approx(-1.0));
}

TEST_CASE("test poly_from_roots complex conjugate") {
    auto roots = std::vector<std::complex<double>>{{0.0, 1.0}, {0.0, -1.0}};
    auto coeffs = poly_from_roots(roots);
    REQUIRE(coeffs.size() == 3);
    CHECK_EQ(coeffs[0], doctest::Approx(1.0));
    CHECK_EQ(coeffs[1], doctest::Approx(0.0));
    CHECK_EQ(coeffs[2], doctest::Approx(1.0));
}

TEST_CASE("test poly_from_roots triple root") {
    auto roots = std::vector<std::complex<double>>{{0.0, 0.0}, {1.0, 0.0}, {-1.0, 0.0}};
    auto coeffs = poly_from_roots(roots);
    REQUIRE(coeffs.size() == 4);
    CHECK_EQ(coeffs[0], doctest::Approx(1.0));
    CHECK_EQ(coeffs[1], doctest::Approx(0.0));
    CHECK_EQ(coeffs[2], doctest::Approx(-1.0));
    CHECK_EQ(coeffs[3], doctest::Approx(0.0));
}

TEST_CASE("test poly_from_roots empty") {
    auto coeffs = poly_from_roots({});
    REQUIRE(coeffs.size() == 1);
    CHECK_EQ(coeffs[0], doctest::Approx(1.0));
}

TEST_CASE("test poly_from_roots aberth reconstruction") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto zs = initial_aberth(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth(h, zs, options);
    REQUIRE(result.second);  // must converge
    auto monic = poly_from_roots(zs);
    REQUIRE(monic.size() == h.size());
    auto scale = h[0];
    for (auto i = 0U; i < h.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(h[i]).epsilon(1e-8));
    }
}

TEST_CASE("test poly_from_roots aberth reconstruction monic") {
    auto h = std::vector<double>{1.0, 0.0, -1.0};
    auto zs = initial_aberth(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth(h, zs, options);
    REQUIRE(result.second);
    auto reconstructed = poly_from_roots(zs);
    REQUIRE(reconstructed.size() == 3);
    CHECK_EQ(reconstructed[0], doctest::Approx(1.0).epsilon(1e-8));
    CHECK_EQ(reconstructed[1], doctest::Approx(0.0).epsilon(1e-8));
    CHECK_EQ(reconstructed[2], doctest::Approx(-1.0).epsilon(1e-8));
}

TEST_CASE("test leja_order empty") {
    auto ordered = leja_order({});
    CHECK(ordered.empty());
}

TEST_CASE("test leja_order single") {
    auto ordered = leja_order({{3.0, 4.0}});
    REQUIRE(ordered.size() == 1);
    CHECK_EQ(ordered[0].real(), doctest::Approx(3.0));
    CHECK_EQ(ordered[0].imag(), doctest::Approx(4.0));
}

TEST_CASE("test leja_order preserves set") {
    using C = std::complex<double>;
    auto points = std::vector<C>{{3.0, 1.0}, {1.0, 2.0}, {0.5, 0.5}, {-2.0, -1.0}};
    auto ordered = leja_order(points);
    REQUIRE(ordered.size() == points.size());
    for (const auto& p : points) {
        auto found = false;
        for (const auto& q : ordered) {
            if (std::abs(p - q) < 1e-14) {
                found = true;
                break;
            }
        }
        CHECK(found);
    }
}

TEST_CASE("test leja_order smallest magnitude first") {
    using C = std::complex<double>;
    auto points = std::vector<C>{{10.0, 0.0}, {0.5, 0.0}, {3.0, 0.0}};
    auto ordered = leja_order(points);
    REQUIRE(ordered.size() == 3);
    CHECK_EQ(ordered[0].real(), doctest::Approx(0.5));  // smallest magnitude first
}

TEST_CASE("test leja_order known 4-point sequence") {
    using C = std::complex<double>;
    auto points = std::vector<C>{{1.0, 0.0}, {-1.0, 0.0}, {0.0, 1.0}, {0.0, -1.0}};
    auto ordered = leja_order(points);
    // All magnitude 1, so any is valid as first. Check set is preserved and size is correct.
    REQUIRE(ordered.size() == 4);
    for (const auto& p : points) {
        auto found = false;
        for (const auto& q : ordered) {
            if (std::abs(p - q) < 1e-14) {
                found = true;
                break;
            }
        }
        CHECK(found);
    }
}

TEST_CASE("test poly_from_autocorr_roots empty") {
    auto coeffs = poly_from_autocorr_roots({});
    REQUIRE(coeffs.size() == 1);
    CHECK_EQ(coeffs[0], doctest::Approx(1.0));
}

TEST_CASE("test poly_from_autocorr_roots reconstruction") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto zs = initial_aberth_autocorr(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth_autocorr(h, zs, options);
    REQUIRE(result.second);
    auto monic = poly_from_autocorr_roots(zs);
    REQUIRE(monic.size() == h.size());
    auto scale = h[0];
    for (auto i = 0U; i < h.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(h[i]).epsilon(1e-8));
    }
}

TEST_CASE("test poly_from_autocorr_roots mt reconstruction") {
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto zs = initial_aberth_autocorr(h);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = aberth_autocorr_mt(h, zs, options);
    REQUIRE(result.second);
    auto monic = poly_from_autocorr_roots(zs);
    REQUIRE(monic.size() == h.size());
    auto scale = h[0];
    for (auto i = 0U; i < h.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(h[i]).epsilon(1e-8));
    }
}

TEST_CASE("test poly_from_autocorr_roots fir") {
    auto r = std::vector<double>{
        -0.00196191, -0.00094597, -0.00023823, 0.00134667,  0.00380494,  0.00681596,  0.0097864,
        0.01186197,  0.0121238,   0.00985211,  0.00474894,  -0.00281751, -0.01173923, -0.0201885,
        -0.02590168, -0.02658216, -0.02035729, -0.00628271, 0.01534627,  0.04279982,  0.0732094,
        0.10275561,  0.12753013,  0.14399228,  0.15265722,  0.14399228,  0.12753013,  0.10275561,
        0.0732094,   0.04279982,  0.01534627,  -0.00628271, -0.02035729, -0.02658216, -0.02590168,
        -0.0201885,  -0.01173923, -0.00281751, 0.00474894,  0.00985211,  0.0121238,   0.01186197,
        0.0097864,   0.00681596,  0.00380494,  0.00134667,  -0.00023823, -0.00094597, -0.00196191};
    auto zs = initial_aberth_autocorr(r);
    auto options = Options();
    options.tolerance = 1e-8;
    auto result = aberth_autocorr(r, zs, options);
    REQUIRE(result.second);
    auto monic = poly_from_autocorr_roots(zs);
    REQUIRE(monic.size() == r.size());
    auto scale = r[0];
    for (auto i = 0U; i < r.size(); ++i) {
        CHECK_EQ(monic[i] * scale, doctest::Approx(r[i]).epsilon(1e-3));
    }
}

TEST_CASE("test leja_order known 5-point sequence") {
    using C = std::complex<double>;
    // Roots with varying magnitudes test the Leja ordering
    auto roots = std::vector<C>{{0.01, 0.0}, {0.1, 0.0}, {1.0, 0.0}, {10.0, 0.0}, {100.0, 0.0}};
    auto ordered = leja_order(roots);
    REQUIRE(ordered.size() == 5);
    // Smallest magnitude first
    CHECK_EQ(ordered[0].real(), doctest::Approx(0.01));
    // All original points present
    for (const auto& r : roots) {
        auto found = false;
        for (const auto& o : ordered) {
            if (std::abs(r - o) < 1e-14) {
                found = true;
                break;
            }
        }
        CHECK(found);
    }
}
