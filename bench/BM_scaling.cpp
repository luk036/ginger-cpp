#include <chrono>
#include <complex>
#include <cstdio>
#include <ginger/aberth.hpp>
#include <ginger/config.hpp>
#include <ginger/rootfinding.hpp>
#include <random>
#include <vector>

using Clock = std::chrono::steady_clock;

static auto real_poly(std::size_t degree, unsigned seed = 42) -> std::vector<double> {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> re(-3.0, 3.0);
    std::uniform_real_distribution<double> im(0.2, 3.0);
    std::vector<std::complex<double>> roots;
    roots.reserve(degree);
    while (roots.size() < degree) {
        const auto a = re(rng);
        const auto b = im(rng);
        roots.emplace_back(a, b);
        if (roots.size() < degree) {
            roots.emplace_back(a, -b);
        }
    }
    return poly_from_roots(roots);
}

static auto max_abs_coeff(const std::vector<double>& c) -> double {
    auto m = 0.0;
    for (const auto v : c) {
        m = std::max(m, std::abs(v));
    }
    return m;
}

static auto scaled(const std::vector<double>& coeffs) -> std::vector<double> {
    const auto scale = ginger::residual_scale(coeffs);
    auto out = coeffs;
    for (auto& c : out) {
        c /= scale;
    }
    return out;
}

static auto poly_residual(const std::vector<double>& coeffs,
                          const std::vector<std::complex<double>>& zs) -> double {
    const auto s = scaled(coeffs);
    auto worst = 0.0;
    for (const auto& z : zs) {
        worst = std::max(worst, std::abs(horner_eval_c(s, z)));
    }
    return worst;
}

static auto factor_residual(const std::vector<double>& coeffs,
                            const std::vector<Vec2>& vrs) -> double {
    const auto s = scaled(coeffs);
    const auto degree = s.size() - 1;
    auto worst = 0.0;
    for (const auto& vr : vrs) {
        auto local = s;
        const auto rem = horner(local, degree, vr);
        worst = std::max(worst, std::max(std::abs(rem.x()), std::abs(rem.y())));
    }
    return worst;
}

template <typename F> static auto best_ms(F&& fn, int repeats = 3) -> double {
    auto best = 1e18;
    for (auto r = 0; r < repeats; ++r) {
        const auto t0 = Clock::now();
        fn();
        const auto dt = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
        best = std::min(best, dt);
    }
    return best;
}

int main() {
    std::printf("%-10s %6s %12s %8s %6s %12s %10s\n", "solver", "degree", "max|coeff|", "niter",
                "conv", "rel.resid", "ms");

    for (const auto degree : {8U, 16U, 32U}) {
        const auto coeffs = real_poly(degree);
        const auto zs = initial_aberth(coeffs);
        const auto ms = best_ms([&] {
            auto local = zs;
            aberth(coeffs, local, ginger::Options{});
        });
        auto local = zs;
        const auto [niter, ok] = aberth(coeffs, local, ginger::Options{});
        std::printf("%-10s %6u %12.3e %8u %6d %12.3e %10.3f\n", "aberth", degree,
                    max_abs_coeff(coeffs), niter, static_cast<int>(ok),
                    poly_residual(coeffs, local), ms);
    }

    for (const auto degree : {8U, 16U, 32U}) {
        const auto coeffs = real_poly(degree);
        const auto vrs = initial_guess(coeffs);
        const auto ms = best_ms([&] {
            auto local = vrs;
            pbairstow_even(coeffs, local, ginger::Options{});
        });
        auto local = vrs;
        const auto [niter, ok] = pbairstow_even(coeffs, local, ginger::Options{});
        std::printf("%-10s %6u %12.3e %8u %6d %12.3e %10.3f\n", "pbairstow", degree,
                    max_abs_coeff(coeffs), niter, static_cast<int>(ok),
                    factor_residual(coeffs, local), ms);
    }

    std::printf("\n%-10s %6s %12s\n", "leja", "n", "ms");
    for (const auto n : {50U, 100U, 200U, 400U}) {
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> d(0.0, 1.0);
        std::vector<std::complex<double>> pts;
        pts.reserve(n);
        for (auto i = 0U; i < n; ++i) {
            pts.emplace_back(d(rng), d(rng));
        }
        const auto ms = best_ms(
            [&] {
                const auto r = leja_order(pts);
                (void)r;
            },
            1);
        std::printf("%-10s %6u %12.3f\n", "leja_order", n, ms);
    }
    return 0;
}
