#pragma once

/**
 * @brief Options
 *
 * The code snippet defines a class called `Options` that represents the options for a specific
 * algorithm or function. It has two public member variables: `max_iters` and `tolerance`.
 */
class Options {
  public:
    unsigned int max_iters = 2000U;
    double tolerance = 1e-14;
};
