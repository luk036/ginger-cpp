/**
 * @file autocorr_mt.hpp
 * @brief Auto-correlation polynomial root-finding (palindromic polynomials) (MT)
 */

#pragma once

#include "autocorr.hpp"

/**
 * @brief Multi-threading Bairstow's method (specific for auto-correlation function)
 *
 * The function `pbairstow_autocorr_mt` implements Bairstow's method for finding the roots of a
 * palindromic (auto-correlation) polynomial using multi-threading.
 *
 * Each thread finds a quadratic factor @f$ x^2 - r_i x - q_i @f$ that respects the palindromic
 * symmetry, where the roots appear in reciprocal pairs:
 * @f[
 *     \begin{bmatrix} \Delta r_i \\ \Delta q_i \end{bmatrix} = -J_i^{-1} \begin{bmatrix} P_i \\ Q_i
 * \end{bmatrix}
 * @f]
 *
 * @dot
 *   digraph pbairstow_ac_flow {
 *     rankdir=LR; bgcolor="transparent";
 *     node [shape=box, style=filled, fillcolor="#d4e6f1"];
 *     coeffs [label="Palindromic\npolynomial", fillcolor="#a9cce3"];
 *     spawn [label="Spawn threads\none per factor", fillcolor="#f9e79f"];
 *     thread [label="Each thread:\nBairstow Newton\n(r_i, q_i) with\nreciprocal pairs",
 * fillcolor="#d4e6f1"]; sync [label="Sync all\nconverged?", shape=diamond, fillcolor="#f9e79f"];
 *     extract [label="Extract\nreciprocal roots\n(r, 1/r)", fillcolor="#d5f5e3"];
 *     done [label="All roots\nfound!", fillcolor="#7fb3d8"];
 *     coeffs -> spawn -> thread -> sync;
 *     sync -> thread [label="No", style=dashed, color="#e74c3c"];
 *     sync -> extract -> done [label="Yes", color="#27ae60"];
 *   }
 * @enddot
 *
 * @param[in] coeffs polynomial
 * @param[in,out] vrs vector of iterates
 * @param[in] options maximum iterations and tolorance
 * @return std::pair<unsigned int, bool>
 */
extern auto pbairstow_autocorr_mt(const std::vector<double>& coeffs, std::vector<Vec2>& vrs,
                                  const ginger::Options& options) -> std::pair<unsigned int, bool>;
