#include <algorithm>                      // for max, min
#include <cstddef>                        // for size_t
#include <ginger/config.hpp>              // for Options
#include <ginger/rootfinding_atomic.hpp>  // for Vec2, pbairstow_even_atomic
#include <utility>                        // for pair
#include <vector>                         // for vector

#include "execution_policy.hpp"  // for ginger::detail::atomic_decoupled_policy, even_bairstow_step

auto pbairstow_even_atomic(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                           const ginger::Options& options) -> std::pair<unsigned int, bool> {
    const auto degree = coeffs.size() - 1;
    ginger::detail::even_bairstow_step step{coeffs, degree, options};
    return ginger::detail::atomic_decoupled_policy::run(vrs, options, step);
}
