#include <benchmark/benchmark.h>
#include <ginger/autocorr.hpp>
#include <ginger/config.hpp>
#include <ginger/rootfinding.hpp>
#include <vector>

auto run_autocorr_st() {
    auto r = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_autocorr(r);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_autocorr_st(r, vrs, options);
    return result;
}

auto run_pbairstow_st() {
    auto r = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_guess(r);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_even_st(r, vrs, options);
    return result;
}

auto run_autocorr_mt() {
    auto r = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_autocorr(r);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_autocorr_mt(r, vrs, options);
    return result;
}

auto run_pbairstow_mt() {
    auto r = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_guess(r);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_even_mt(r, vrs, options);
    return result;
}

static void Autocorr_ST(benchmark::State& state) {
    for (auto _ : state) run_autocorr_st();
}
BENCHMARK(Autocorr_ST);

static void PBairstow_ST(benchmark::State& state) {
    for (auto _ : state) run_pbairstow_st();
}
BENCHMARK(PBairstow_ST);

static void Autocorr_MT(benchmark::State& state) {
    for (auto _ : state) run_autocorr_mt();
}
BENCHMARK(Autocorr_MT);

static void PBairstow_MT(benchmark::State& state) {
    for (auto _ : state) run_pbairstow_mt();
}
BENCHMARK(PBairstow_MT);

BENCHMARK_MAIN();
