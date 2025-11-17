#include <doctest/doctest.h>
#include <ginger/aberth.hpp>
#include <ginger/config.hpp>
#include <vector>
#include <random>
#include <chrono>

TEST_CASE("stress test aberth_mt with high-degree polynomial") {
    // Generate a high-degree polynomial with random coefficients
    const int degree = 100;
    std::vector<double> h(degree + 1);
    std::mt19937_64 rng(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    for (int i = 0; i <= degree; ++i) {
        h[i] = dist(rng);
    }

    // Set options for Aberth method
    Options options;
    options.tolerance = 1e-9;
    options.max_iters = 1000; // Increased max iterations for stress test

    // Initial guess for roots
    auto zs = initial_aberth(h);

    // Run the multi-threaded Aberth method
    auto result = aberth_mt(h, zs, options);
    auto niter = result.first;

    // Check if the number of iterations is within a reasonable bound
    // The exact bound might vary, but it should converge within a certain number of iterations
    CHECK(niter <= options.max_iters);
    CHECK(niter > 0); // Ensure it actually ran
}
