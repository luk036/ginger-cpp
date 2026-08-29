/**
 * @file rootfinding_atomic.hpp
 * @brief Parallel Bairstow root-finding methods for real polynomials (atomic)
 */

#pragma once

#include "rootfinding.hpp"

/**
 * @brief Atomic multi-threading Bairstow's method (even degree only)
 *
 * The `pbairstow_even_atomic` function implements Bairstow's method for finding the roots of a
 * real polynomial with an even degree using multi-threading with an atomic working buffer.
 *
 * Unlike the snapshot-based Jacobi variant (`pbairstow_even_mt`), a single atomic working buffer
 * is built once and reused across all iterations: each thread owns exactly one slot (single-writer,
 * multi-reader) and reads the latest values of the other slots while iterating. The whole
 * @f$ (r_i, q_i) @f$ pair is read and written atomically via `std::atomic<Vec2>` (16 bytes, not
 * lock-free on MSVC but correct) so that `suppress_old`/`delta_scalar` always consume one coherent
 * pair. Iterations become asynchronous/in-place (Gauss-Seidel-like): each thread runs its own
 * iteration loop independently with no per-iteration synchronization and exits when its own factors
 * converge or the maximum number of iterations is exceeded. The iteration count is therefore
 * NON-DETERMINISTIC (the returned count is the maximum across threads) — convergence, not a fixed
 * iteration count, must be asserted.
 *
 * Each thread handles one quadratic factor @f$ x^2 - r_i x - q_i @f$, applying:
 * @f[
 *     \begin{bmatrix} \Delta r_i \\ \Delta q_i \end{bmatrix} = -J_i^{-1} \begin{bmatrix} P_i \\ Q_i
 * \end{bmatrix}
 * @f]
 * where @f$ P_i, Q_i @f$ are the remainders from synthetic division.
 *
 * @dot
 *   digraph pbairstow_atomic_flow {
 *     bgcolor="transparent";
 *     rankdir=LR;
 *     node [shape=box, style=filled, fillcolor="#d4e6f1"];
 *     init [label="Initial\nquadratic factors", fillcolor="#a9cce3"];
 *     buffer [label="Atomic buffer\nbuilt once", fillcolor="#f9e79f"];
 *     spawn [label="Spawn threads\none per factor"];
 *     thread [label="Each thread:\nread latest slots,\nwrite own slot", fillcolor="#d4e6f1"];
 *     sync [label="Sync +\ncheck all\nconverged?", shape=diamond, fillcolor="#f9e79f"];
 *     extract [label="Extract\nfactors", fillcolor="#d5f5e3"];
 *     done [label="All roots\nfound!", fillcolor="#7fb3d8"];
 *     init -> buffer -> spawn -> thread -> sync;
 *     sync -> thread [label="No", style=dashed, color="#e74c3c"];
 *     sync -> extract -> done [label="Yes", color="#27ae60"];
 *   }
 * @enddot
 *
 * @param[in] coeffs The `coeffs` parameter is a vector representing the coefficients of the
 * polynomial. Each element of the vector corresponds to the coefficient of a term in the
 * polynomial, starting from the highest degree term and ending with the constant term. For example,
 * if the polynomial is `3x^2 + 2x + 1`, the coefficients vector would be `{3, 2, 1}`.
 * @param[in,out] vrs `vrs` is a vector of iterates, which represents the initial guesses for the
 * roots of the polynomial. Bairstow's method will update these iterates iteratively until the
 * desired tolerance is reached or the maximum number of iterations is reached.
 * @param[in] options The `options` parameter is an object of type `Options` which contains the
 * maximum number of iterations (`max_iters`) and the tolerance (`tolerance`). These options are
 * used to control the convergence criteria for Bairstow's method.
 *
 * @return The function `pbairstow_even_atomic` returns a `std::pair<unsigned int, bool>`. The first
 * element of the pair represents the number of iterations performed, and the second element
 * represents whether the method converged to a solution within the specified tolerance.
 */
extern auto pbairstow_even_atomic(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                                  const ginger::Options& options) -> std::pair<unsigned int, bool>;
