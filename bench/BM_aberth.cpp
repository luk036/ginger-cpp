#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include <ginger/aberth.hpp>
#include <ginger/aberth_mt.hpp>
#include <ginger/config.hpp>
#include <vector>

static const auto global_r = std::vector<double>{
    -0.00196191, -0.00094597, -0.00023823, 0.00134667,  0.00380494,  0.00681596,  0.0097864,
    0.01186197,  0.0121238,   0.00985211,  0.00474894,  -0.00281751, -0.01173923, -0.0201885,
    -0.02590168, -0.02658216, -0.02035729, -0.00628271, 0.01534627,  0.04279982,  0.0732094,
    0.10275561,  0.12753013,  0.14399228,  0.15265722,  0.14399228,  0.12753013,  0.10275561,
    0.0732094,   0.04279982,  0.01534627,  -0.00628271, -0.02035729, -0.02658216, -0.02590168,
    -0.0201885,  -0.01173923, -0.00281751, 0.00474894,  0.00985211,  0.0121238,   0.01186197,
    0.0097864,   0.00681596,  0.00380494,  0.00134667,  -0.00023823, -0.00094597, -0.00196191};

static const auto degree8
    = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};

int main() {
    {
        ankerl::nanobench::Bench bench;
        bench.title("Aberth root-finding (ST)").unit("op").warmup(100).epochs(50).minEpochIterations(20);

        bench.run("FIR_Aberth", [&] {
            auto vrs = initial_aberth(global_r);
            Options opts;
            opts.tolerance = 1e-8;
            auto result = aberth(global_r, vrs, opts);
            ankerl::nanobench::doNotOptimizeAway(result);
        });

        bench.run("Aberth_ST", [&] {
            auto vrs = initial_aberth(degree8);
            Options opts;
            opts.tolerance = 1e-12;
            auto result = aberth(degree8, vrs, opts);
            ankerl::nanobench::doNotOptimizeAway(result);
        });
    }

    {
        ankerl::nanobench::Bench bench;
        bench.title("Aberth root-finding (MT)")
            .unit("op")
            .warmup(10)
            .epochs(30)
            .minEpochIterations(400);

        bench.run("FIR_Aberth_MT", [&] {
            auto vrs = initial_aberth(global_r);
            Options opts;
            opts.tolerance = 1e-8;
            auto result = aberth_mt(global_r, vrs, opts);
            ankerl::nanobench::doNotOptimizeAway(result);
        });

        bench.run("Aberth_MT", [&] {
            auto vrs = initial_aberth(degree8);
            Options opts;
            opts.tolerance = 1e-12;
            auto result = aberth_mt(degree8, vrs, opts);
            ankerl::nanobench::doNotOptimizeAway(result);
        });
    }
}
