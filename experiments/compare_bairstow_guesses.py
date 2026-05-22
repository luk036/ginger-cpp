"""
Experiments comparing LDS vs naive initial guess strategies in Bairstow's method.

Claims to verify:
1. LDS and naive generate different angle sets for Bairstow guess generators
2. LDS converges in fewer iterations on the ST path (Gauss-Seidel)
3. Ordering WITHIN the vector matters for ST path
"""
import math
import random
import statistics
import sys
from typing import Callable, List, Tuple

from lds_gen.lds import VdCorput
from mywheel.robin import Robin

from ginger.rootfinding import (
    Options,
    delta,
    horner,
    initial_guess,
    suppress,
)
from ginger.autocorr import pbairstow_autocorr
from ginger.vector2 import Vector2

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8")


# ─── Naive (non-LDS) versions for comparison ─────────────────────────────────

def initial_guess_naive(coeffs: List[float]) -> List[Vector2]:
    """Naive equidistant version of initial_guess."""
    degree = len(coeffs) - 1
    center = -coeffs[1] / (degree * coeffs[0])
    poly_c = sum(c * (center ** (degree - i)) for i, c in enumerate(coeffs))
    radius = pow(abs(poly_c), 1.0 / degree)
    quad_term = center * center + radius * radius
    deg = degree // 2 * 2
    k = math.pi / deg
    num = deg // 2
    return [
        Vector2(2 * (center + radius * math.cos(k * (2 * i + 1))),
                -(quad_term + 2 * center * radius * math.cos(k * (2 * i + 1))))
        for i in range(num)
    ]


def initial_autocorr_naive(coeffs: List[float]) -> List[Vector2]:
    """Naive equidistant version of initial_autocorr."""
    degree = len(coeffs) - 1
    radius = pow(abs(coeffs[-1]), 1.0 / degree)
    if radius < 1:
        radius = 1 / radius
    degree //= 2
    angle_step = math.pi / degree
    quad_term = radius * radius
    return [
        Vector2(2 * radius * math.cos(angle_step * i), -quad_term)
        for i in range(1, degree, 2)
    ]


def initial_autocorr_lds(coeffs: List[float]) -> List[Vector2]:
    """LDS version of initial_autocorr."""
    degree = len(coeffs) - 1
    radius = pow(abs(coeffs[-1]), 1.0 / degree)
    if radius < 1:
        radius = 1 / radius
    degree //= 2
    num_points = degree // 2
    quad_term = radius * radius
    vgen = VdCorput(2)
    vgen.reseed(1)
    return [
        Vector2(2 * radius * math.cos(math.pi * vgen.pop()), -quad_term)
        for _ in range(num_points)
    ]


# ─── Test polynomials ────────────────────────────────────────────────────────

P1 = [10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0]   # deg 8 (palindromic)

P2 = [5.0, 2.0, 9.0, 6.0, 2.0]                                     # deg 4

# FIR filter, degree 48 (from test_aberth.py)
P3 = [
    -0.00196191, -0.00094597, -0.00023823, 0.00134667,  0.00380494,  0.00681596,
     0.0097864,   0.01186197,  0.0121238,   0.00985211,  0.00474894,  -0.00281751,
    -0.01173923, -0.0201885,  -0.02590168, -0.02658216, -0.02035729, -0.00628271,
     0.01534627,  0.04279982,  0.0732094,   0.10275561,  0.12753013,  0.14399228,
     0.15265722,  0.14399228,  0.12753013,  0.10275561,  0.0732094,   0.04279982,
     0.01534627, -0.00628271, -0.02035729, -0.02658216, -0.02590168, -0.0201885,
    -0.01173923, -0.00281751,  0.00474894,  0.00985211,  0.0121238,   0.01186197,
     0.0097864,   0.00681596,  0.00380494,  0.00134667, -0.00023823, -0.00094597,
    -0.00196191,
]

# Synthetic degree 20: alternating coefficients with 3x growth factor
P4 = [1.0]
for k in range(20):
    P4.append(P4[-1] * (1.5 if k % 2 == 0 else 0.5))

ALL_POLYS = [
    ("P2 (deg 4)", P2),
    ("P1 (deg 8)", P1),
    ("P4 (deg 20)", P4),
    ("P3 (deg 48 FIR)", P3),
]

ALL_POLYS_BIG = ALL_POLYS


# ─── Angle distribution ─────────────────────────────────────────────────────

def angles_from_guess_naive(n: int) -> List[float]:
    """Angles for initial_guess with naive distribution."""
    deg = n
    k = math.pi / deg
    num = deg // 2
    return [k * (2 * i + 1) for i in range(num)]


def angles_from_guess_lds(n: int) -> List[float]:
    """Angles for initial_guess with LDS distribution."""
    vgen = VdCorput(2)
    vgen.reseed(1)
    deg = n
    num = deg // 2
    return [math.pi * vgen.pop() for _ in range(num)]


def angles_from_autocorr_naive(n: int) -> List[float]:
    """Angles for initial_autocorr with naive distribution."""
    deg = n // 2
    angle_step = math.pi / deg
    return [angle_step * i for i in range(1, deg, 2)]


def angles_from_autocorr_lds(n: int) -> List[float]:
    """Angles for initial_autocorr with LDS distribution."""
    vgen = VdCorput(2)
    vgen.reseed(1)
    deg = n // 2
    num = deg // 2
    return [math.pi * vgen.pop() for _ in range(num)]


def exp1_angle_distribution() -> None:
    """Show LDS and naive produce DIFFERENT angle sets for Bairstow."""
    print("=" * 72)
    print("Experiment 1: Bairstow initial guess angle distribution")
    print("=" * 72)

    for n in [8, 12, 16]:
        g_naive = angles_from_guess_naive(n)
        g_lds = angles_from_guess_lds(n)
        a_naive = angles_from_autocorr_naive(n)
        a_lds = angles_from_autocorr_lds(n)

        g_same = all(abs(a - b) < 1e-10 for a, b in zip(
            sorted(g_naive), sorted(g_lds)))
        a_same = all(abs(a - b) < 1e-10 for a, b in zip(
            sorted(a_naive), sorted(a_lds)))

        print(f"\ndegree={n}:")
        print(f"  initial_guess naive: {[f'{a*180/math.pi:.0f}' for a in sorted(g_naive)]}")
        print(f"  initial_guess LDS:   {[f'{a*180/math.pi:.0f}' for a in sorted(g_lds)]}")
        print(f"  Sets identical: {g_same}")
        print(f"  initial_autocorr naive: {[f'{a*180/math.pi:.0f}' for a in sorted(a_naive)]}")
        print(f"  initial_autocorr LDS:   {[f'{a*180/math.pi:.0f}' for a in sorted(a_lds)]}")
        print(f"  Sets identical: {a_same}")


# ─── Tolerance helper ────────────────────────────────────────────────────────

def get_tol(coeffs: List[float]) -> float:
    deg = len(coeffs) - 1
    if deg >= 40:  return 1e-8
    if deg >= 20:  return 1e-10
    return 1e-12


# ─── Bairstow ST convergence ─────────────────────────────────────────────────

def run_pbairstow(
    coeffs: List[float], vrs: List[Vector2], tol: float | None = None, max_iters: int = 2000
) -> Tuple[List[Vector2], int, bool]:
    """Simplified ST Bairstow (single-thread, sequential updates)."""
    if tol is None:
        tol = get_tol(coeffs)
    num_factors = len(vrs)
    degree = len(coeffs) - 1
    robin = Robin(num_factors)
    for niter in range(max_iters):
        tolerance = 0.0
        for i, vri in enumerate(vrs):
            coeffs1 = coeffs.copy()
            vA = horner(coeffs1, degree, vri)
            tol_i = max(abs(vA.x), abs(vA.y))
            vA1 = horner(coeffs1, degree - 2, vri)
            tolerance = max(tol_i, tolerance)
            for j in robin.exclude(i):
                vA, vA1 = suppress(vA, vA1, vri, vrs[j])
            vrs[i] -= delta(vA, vri, vA1)
        if tolerance < tol:
            return vrs, niter, True
    return vrs, max_iters, False


def exp2_pbairstow_convergence() -> None:
    """Compare Bairstow ST convergence: LDS vs naive."""
    print("\n" + "=" * 72)
    print("Experiment 2: Bairstow ST convergence (LDS vs naive)")
    print("=" * 72)
    for pname, coeffs in ALL_POLYS:
        vrs_lds = initial_guess(coeffs)
        _, n_lds, f_lds = run_pbairstow(coeffs, vrs_lds[:])
        vrs_naive = initial_guess_naive(coeffs)
        _, n_naive, f_naive = run_pbairstow(coeffs, vrs_naive[:])
        print(f"\n{pname}:")
        print(f"  LDS:   {n_lds:4d} iters [{ 'OK' if f_lds else 'FAIL' }]")
        print(f"  Naive: {n_naive:4d} iters [{ 'OK' if f_naive else 'FAIL' }]")
        if f_lds and f_naive and n_naive > 0:
            diff = n_naive - n_lds
            if diff > 0:
                print(f"  -> LDS {diff} iters FASTER ({diff/n_naive*100:.0f}%)")
            elif diff == 0:
                print(f"  -> Same iterations")
            else:
                print(f"  -> Naive {-diff} iters FASTER")


# ─── Autocorr Bairstow convergence ───────────────────────────────────────────

def run_pbairstow_autocorr(
    coeffs: List[float], vrs: List[Vector2], tol: float | None = None, max_iters: int = 2000
) -> Tuple[List[Vector2], int, bool]:
    """Simplified ST autocorr Bairstow."""
    if tol is None:
        tol = get_tol(coeffs)
    num_factors = len(vrs)
    degree = len(coeffs) - 1
    robin = Robin(num_factors)
    for niter in range(max_iters):
        tolerance = 0.0
        for i, vri in enumerate(vrs):
            coeffs1 = coeffs.copy()
            vA = horner(coeffs1, degree, vri)
            tol_i = max(abs(vA.x), abs(vA.y))
            vA1 = horner(coeffs1, degree - 2, vri)
            tolerance = max(tol_i, tolerance)
            for j in robin.exclude(i):
                vrj = vrs[j]
                vA, vA1 = suppress(vA, vA1, vri, vrj)
                vrn = Vector2(-vrj.x, 1.0) / vrj.y
                vA, vA1 = suppress(vA, vA1, vri, vrn)
            vrin = Vector2(-vri.x, 1.0) / vri.y
            vA, vA1 = suppress(vA, vA1, vri, vrin)
            vrs[i] -= delta(vA, vri, vA1)
        if tolerance < tol:
            return vrs, niter, True
    return vrs, max_iters, False


def exp3_autocorr_convergence() -> None:
    """Compare autocorr Bairstow convergence: LDS vs naive."""
    print("\n" + "=" * 72)
    print("Experiment 3: Autocorr Bairstow convergence (LDS vs naive)")
    print("=" * 72)
    for pname, coeffs in ALL_POLYS:
        if len(coeffs) < 4:
            print(f"\n{pname}: skip (need even degree >= 4)")
            continue
        vrs_l = initial_autocorr_lds(coeffs)
        _, n_l, f_l = run_pbairstow_autocorr(coeffs, vrs_l[:])
        vrs_n = initial_autocorr_naive(coeffs)
        _, n_n, f_n = run_pbairstow_autocorr(coeffs, vrs_n[:])
        print(f"\n{pname}:")
        print(f"  LDS:   {n_l:4d} iters [{ 'OK' if f_l else 'FAIL' }]")
        print(f"  Naive: {n_n:4d} iters [{ 'OK' if f_n else 'FAIL' }]")
        if f_l and f_n and n_n > 0:
            diff = n_n - n_l
            if diff > 0:
                print(f"  -> LDS {diff} iters FASTER ({diff/n_n*100:.0f}%)")
            elif diff == 0:
                print(f"  -> Same iterations")
            else:
                print(f"  -> Naive {-diff} iters FASTER")


# ─── Ordering effect (shuffle test) ──────────────────────────────────────────

def exp4_ordering_effect() -> None:
    """Test whether order matters for Bairstow and autocorr."""
    print("\n" + "=" * 72)
    print("Experiment 4: Ordering effect (shuffle test)")
    print("=" * 72)

    for label, coeffs, guess_fn, runner_fn in [
        ("Bairstow", P2, initial_guess, run_pbairstow),
        ("Bairstow", P1, initial_guess, run_pbairstow),
        ("Bairstow", P4, initial_guess, run_pbairstow),
    ]:
        _, n_orig, f_orig = runner_fn(coeffs, guess_fn(coeffs))
        pname = f"{label} (deg {len(coeffs)-1})"
        if not f_orig:
            print(f"\n{pname}: did not converge, skipping")
            continue
        random.seed(42)
        shuffle_results = []
        for _ in range(10):
            shuffled = guess_fn(coeffs)
            random.shuffle(shuffled)
            _, n, f = runner_fn(coeffs, shuffled)
            shuffle_results.append((n, f))
        converged = [n for n, f in shuffle_results if f]
        print(f"\n{pname}:")
        print(f"  Original LDS order: {n_orig} iters")
        if converged:
            print(f"  Shuffled (10 trials): min={min(converged)}, "
                  f"max={max(converged)}, mean={statistics.mean(converged):.1f}")
            worst = max(converged) - n_orig
            pct = worst / max(converged) * 100 if max(converged) > 0 else 0
            print(f"  -> Worst shuffle slowdown: {worst} iters ({pct:.0f}%)")

    for label, coeffs, guess_fn, runner_fn in [
        ("Autocorr", P1, initial_autocorr_lds, run_pbairstow_autocorr),
        ("Autocorr", P4, initial_autocorr_lds, run_pbairstow_autocorr),
    ]:
        vrs_ref = guess_fn(coeffs)
        _, n_orig, f_orig = runner_fn(coeffs, guess_fn(coeffs))
        if not f_orig:
            print(f"\n{label}: did not converge, skipping")
            continue
        random.seed(42)
        shuffle_results = []
        for _ in range(10):
            shuffled = guess_fn(coeffs)
            random.shuffle(shuffled)
            _, n, f = runner_fn(coeffs, shuffled)
            shuffle_results.append((n, f))
        converged = [n for n, f in shuffle_results if f]
        print(f"\n{label}:")
        print(f"  Original LDS order: {n_orig} iters")
        if converged:
            print(f"  Shuffled (10 trials): min={min(converged)}, "
                  f"max={max(converged)}, mean={statistics.mean(converged):.1f}")
            worst = max(converged) - n_orig
            pct = worst / max(converged) * 100 if max(converged) > 0 else 0
            print(f"  -> Worst shuffle slowdown: {worst} iters ({pct:.0f}%)")


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("Bairstow Initial Guess Strategy Experiments")
    print("=" * 72)
    exp1_angle_distribution()
    exp2_pbairstow_convergence()
    exp3_autocorr_convergence()
    exp4_ordering_effect()
    print("\n" + "=" * 72)
    print("Done.")


if __name__ == "__main__":
    main()
