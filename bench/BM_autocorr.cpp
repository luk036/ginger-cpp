#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

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

int main() {
    {
        ankerl::nanobench::Bench bench;
        bench.title("Autocorr root-finding")
            .unit("op")
            .warmup(100)
            .epochs(50)
            .minEpochIterations(50000);

        bench.run("Autocorr_ST", [&] {
            auto result = run_autocorr_st();
            ankerl::nanobench::doNotOptimizeAway(result);
        });

        bench.run("PBairstow_ST", [&] {
            auto result = run_pbairstow_st();
            ankerl::nanobench::doNotOptimizeAway(result);
        });

        bench.run("Autocorr_MT", [&] {
            auto result = run_autocorr_mt();
            ankerl::nanobench::doNotOptimizeAway(result);
        });
    }

    // PBairstow_MT is slower (~360µs/iter), so use lower minEpochIterations
    {
        ankerl::nanobench::Bench bench;
        bench.title("Autocorr root-finding (PBairstow_MT)")
            .unit("op")
            .warmup(10)
            .epochs(30)
            .minEpochIterations(2000);

        bench.run("PBairstow_MT", [&] {
            auto result = run_pbairstow_mt();
            ankerl::nanobench::doNotOptimizeAway(result);
        });
    }
}
