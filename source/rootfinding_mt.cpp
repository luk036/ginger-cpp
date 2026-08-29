#include <algorithm>                  // for max, min
#include <cstddef>                    // for size_t
#include <ginger/config.hpp>          // for Options
#include <ginger/rootfinding_mt.hpp>  // for Vec2, pbairstow_even_mt
#include <utility>                    // for pair
#include <vector>                     // for vector

#include "execution_policy.hpp"  // for ginger::detail::even_bairstow_step, jacobi_mt_policy

auto pbairstow_even_mt(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                       const ginger::Options& options) -> std::pair<unsigned int, bool> {
    const auto degree = coeffs.size() - 1;
    ginger::detail::even_bairstow_step step{coeffs, degree, options};
    return ginger::detail::jacobi_mt_policy::run(vrs, options, step);
}
