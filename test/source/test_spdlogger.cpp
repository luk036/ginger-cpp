#include <doctest/doctest.h>
#include "ginger/logger.hpp"
#include "ginger/aberth.hpp"
#include "ginger/config.hpp"

TEST_CASE("Spdlogger integration test") {
    // Test wrapper function
    ginger::log_with_spdlog("Test message 1");
    ginger::log_with_spdlog("Test message 2");
    ginger::log_with_spdlog("Test message 3");
    
    // Test with polynomial root finding
    std::vector<double> coeffs = {1.0, -6.0, 11.0, -6.0}; // (x-1)(x-2)(x-3)
    auto initial = initial_aberth(coeffs);
    Options options;
    
    ginger::log_with_spdlog("Starting polynomial root finding");
    auto [iters, converged] = aberth(coeffs, initial, options);
    
    if (converged) {
        ginger::log_with_spdlog("Polynomial root finding converged");
    } else {
        ginger::log_with_spdlog("Polynomial root finding did not converge");
    }
    
    CHECK(true); // Test passes if we get here without crashes
}