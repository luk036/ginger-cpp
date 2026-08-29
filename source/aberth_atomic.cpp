#include <algorithm>                 // for max, min
#include <cstddef>                   // for size_t
#include <ginger/aberth_atomic.hpp>  // for aberth_atomic, aberth_autocorr_atomic
#include <ginger/config.hpp>         // for Options
#include <utility>                   // for pair
#include <vector>                    // for vector

#include "execution_policy.hpp"  // for ginger::detail::aberth_autocorr_step, aberth_step

using std::vector;
using Complex = std::complex<double>;

auto aberth_atomic(const vector<double>& coeffs, vector<Complex>& zs,
                   const ginger::Options& options = ginger::Options())
    -> std::pair<unsigned int, bool> {
    auto coeffs1 = ginger::detail::derivative_coeffs(coeffs);
    ginger::detail::aberth_step step{coeffs, coeffs1};
    return ginger::detail::atomic_decoupled_policy::run(zs, options, step);
}

auto aberth_autocorr_atomic(const vector<double>& coeffs, vector<Complex>& zs,
                            const ginger::Options& options = ginger::Options())
    -> std::pair<unsigned int, bool> {
    auto coeffs1 = ginger::detail::derivative_coeffs(coeffs);
    ginger::detail::aberth_autocorr_step step{coeffs, coeffs1};
    return ginger::detail::atomic_decoupled_policy::run(zs, options, step);
}
