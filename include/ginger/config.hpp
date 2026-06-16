/**
 * @file config.hpp
 * @brief Algorithm options (max iterations, tolerance)
 */

#pragma once

/**
 * @brief Options for convergence-based algorithms
 *
 * Stores maximum iteration count and tolerance for convergence checks
 * used by Bairstow and Aberth root-finding methods.
 */
class Options {
  public:
    unsigned int max_iters = 2000U;
    double tolerance = 1e-14;
};
