#include <ginger/aberth.hpp>
#include <ginger/config.hpp>
#include <ginger/version.h>

#include <cxxopts.hpp>
#include <iostream>
#include <string>
#include <vector>

auto main(int argc, char** argv) -> int {
    cxxopts::Options options("Ginger", "Polynomial root-finding demo");
    options.add_options()("h,help", "Print usage")("v,version", "Print version");

    const auto result = options.parse(argc, argv);
    if (result.count("help") > 0) {
        std::cout << options.help() << '\n';
        return 0;
    }
    if (result.count("version") > 0) {
        std::cout << "Ginger, version " << GINGER_VERSION << '\n';
        return 0;
    }

    const std::vector<double> coeffs{1.0, -3.0, 2.0};
    auto zs = initial_aberth(coeffs);
    const auto [iters, converged] = aberth(coeffs, zs, Options{});

    std::cout << "Ginger: roots of x^2 - 3x + 2 (converged=" << std::boolalpha << converged
              << ", iters=" << iters << "):";
    for (const auto& z : zs) {
        std::cout << ' ' << z;
    }
    std::cout << '\n';

    return 0;
}
