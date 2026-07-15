/**
 * @file config.hpp
 * @brief Algorithm options (max iterations, tolerance)
 */

#pragma once

/**
 * @brief Options for convergence-based algorithms
 *
 * Stores maximum iteration count, tolerance for global convergence,
 * and per-root tolerance for convergence checks
 * used by Bairstow and Aberth root-finding methods.
 */
class Options {
  public:
    unsigned int max_iters = 2000U;
    double tolerance = 1e-12;
    double tol_ind = 1e-15;
};
