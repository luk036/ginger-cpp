#include <algorithm>               // for max, min
#include <cstddef>                 // for size_t
#include <ginger/autocorr_mt.hpp>  // for Vec2, pbairstow_autocorr_mt
#include <ginger/config.hpp>       // for Options
#include <utility>                 // for pair
#include <vector>                  // for vector

#include "execution_policy.hpp"  // for ginger::detail::autocorr_bairstow_step, jacobi_mt_policy

auto pbairstow_autocorr_mt(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                           const ginger::Options& options) -> std::pair<unsigned int, bool> {
    const auto degree = coeffs.size() - 1;
    ginger::detail::autocorr_bairstow_step step{coeffs, degree, options};
    return ginger::detail::jacobi_mt_policy::run(vrs, options, step);
}
