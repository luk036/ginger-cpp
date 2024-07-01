#include <benchmark/benchmark.h>  // for BENCHMARK, State, BENCHMARK_MAIN

#include <ginger/autocorr.hpp>     // for initial_autocorr, pbairstow_auto...
#include <ginger/config.hpp>       // for Options
#include <ginger/rootfinding.hpp>  // for Options, initial_guess, pbairsto...
#include <vector>                  // for vector

static const auto global_r = std::vector<double>{
    -0.00196191, -0.00094597, -0.00023823, 0.00134667,  0.00380494,  0.00681596,  0.0097864,
    0.01186197,  0.0121238,   0.00985211,  0.00474894,  -0.00281751, -0.01173923, -0.0201885,
    -0.02590168, -0.02658216, -0.02035729, -0.00628271, 0.01534627,  0.04279982,  0.0732094,
    0.10275561,  0.12753013,  0.14399228,  0.15265722,  0.14399228,  0.12753013,  0.10275561,
    0.0732094,   0.04279982,  0.01534627,  -0.00628271, -0.02035729, -0.02658216, -0.02590168,
    -0.0201885,  -0.01173923, -0.00281751, 0.00474894,  0.00985211,  0.0121238,   0.01186197,
    0.0097864,   0.00681596,  0.00380494,  0.00134667,  -0.00023823, -0.00094597, -0.00196191};

auto run_fir_autocorr() -> std::pair<unsigned int, bool> {
    auto r = global_r;
    auto vrs = initial_autocorr(r);
    auto options = Options();
    options.tolerance = 1e-2;
    auto result = pbairstow_autocorr(r, vrs, options);
    return result;
}

auto run_fir_pbairstow() -> std::pair<unsigned int, bool> {
    // auto r = std::vector<double>{
    //     -0.00196191, -0.00094597, -0.00023823, 0.00134667,  0.00380494,  0.00681596,  0.0097864,
    //     0.01186197,  0.0121238,   0.00985211,  0.00474894,  -0.00281751, -0.01173923, -0.0201885,
    //     -0.02590168, -0.02658216, -0.02035729, -0.00628271, 0.01534627,  0.04279982,  0.0732094,
    //     0.10275561,  0.12753013,  0.14399228,  0.15265722,  0.14399228,  0.12753013,  0.10275561,
    //     0.0732094,   0.04279982,  0.01534627,  -0.00628271, -0.02035729, -0.02658216,
    //     -0.02590168, -0.0201885,  -0.01173923, -0.00281751, 0.00474894,  0.00985211,  0.0121238,
    //     0.01186197, 0.0097864,   0.00681596,  0.00380494,  0.00134667,  -0.00023823, -0.00094597,
    //     -0.00196191};
    auto r = global_r;
    auto vrs = initial_guess(r);
    auto options = Options();
    options.tolerance = 1e-2;
    auto result = pbairstow_even(r, vrs, options);
    return result;
}

/**
 * @brief
 *
 * @param[in,out] state
 */
static void FIR_Autocorr(benchmark::State &state) {
    while (state.KeepRunning()) {
        auto result = run_fir_autocorr();
        benchmark::DoNotOptimize(result);
    }
}

// Register the function as a benchmark
BENCHMARK(FIR_Autocorr);

/**
 * @brief
 *
 * @param[in,out] state
 */
static void FIR_PBairstow(benchmark::State &state) {
    while (state.KeepRunning()) {
        auto result = run_fir_pbairstow();
        benchmark::DoNotOptimize(result);
    }
}

// Register the function as a benchmark
BENCHMARK(FIR_PBairstow);

BENCHMARK_MAIN();
