#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include <ginger/autocorr.hpp>
#include <ginger/autocorr_atomic.hpp>
#include <ginger/autocorr_mt.hpp>
#include <ginger/config.hpp>
#include <ginger/rootfinding.hpp>
#include <ginger/rootfinding_mt.hpp>
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

auto run_autocorr_atomic() {
    auto r = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_autocorr(r);
    auto options = Options();
    options.tolerance = 1e-12;
    auto result = pbairstow_autocorr_atomic(r, vrs, options);
    return result;
}

static const auto global_r = std::vector<double>{
    -0.00196191, -0.00094597, -0.00023823, 0.00134667,  0.00380494,  0.00681596,  0.0097864,
    0.01186197,  0.0121238,   0.00985211,  0.00474894,  -0.00281751, -0.01173923, -0.0201885,
    -0.02590168, -0.02658216, -0.02035729, -0.00628271, 0.01534627,  0.04279982,  0.0732094,
    0.10275561,  0.12753013,  0.14399228,  0.15265722,  0.14399228,  0.12753013,  0.10275561,
    0.0732094,   0.04279982,  0.01534627,  -0.00628271, -0.02035729, -0.02658216, -0.02590168,
    -0.0201885,  -0.01173923, -0.00281751, 0.00474894,  0.00985211,  0.0121238,   0.01186197,
    0.0097864,   0.00681596,  0.00380494,  0.00134667,  -0.00023823, -0.00094597, -0.00196191};

auto run_fir_autocorr_mt() {
    auto r = global_r;
    auto vrs = initial_autocorr(r);
    auto options = Options();
    options.tolerance = 1e-2;
    auto result = pbairstow_autocorr_mt(r, vrs, options);
    return result;
}

auto run_fir_autocorr_atomic() {
    auto r = global_r;
    auto vrs = initial_autocorr(r);
    auto options = Options();
    options.tolerance = 1e-2;
    auto result = pbairstow_autocorr_atomic(r, vrs, options);
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

    // PBairstow_MT and Autocorr_Atomic are slower on d8 (dispatch / MSVC atomic
    // spinlock), so use lower minEpochIterations
    {
        ankerl::nanobench::Bench bench;
        bench.title("Autocorr root-finding (PBairstow_MT / Autocorr_Atomic)")
            .unit("op")
            .warmup(10)
            .epochs(30)
            .minEpochIterations(2000);

        bench.run("PBairstow_MT", [&] {
            auto result = run_pbairstow_mt();
            ankerl::nanobench::doNotOptimizeAway(result);
        });

        bench.run("Autocorr_Atomic", [&] {
            auto result = run_autocorr_atomic();
            ankerl::nanobench::doNotOptimizeAway(result);
        });
    }

    // FIR (48-tap, 12 factors): the real multi-threaded autocorr path
    {
        ankerl::nanobench::Bench bench;
        bench.title("Autocorr root-finding FIR (MT vs Atomic)")
            .unit("op")
            .warmup(10)
            .epochs(30)
            .minEpochIterations(200);

        bench.run("FIR_AutoCorr_MT", [&] {
            auto result = run_fir_autocorr_mt();
            ankerl::nanobench::doNotOptimizeAway(result);
        });

        bench.run("FIR_AutoCorr_Atomic", [&] {
            auto result = run_fir_autocorr_atomic();
            ankerl::nanobench::doNotOptimizeAway(result);
        });
    }
}
