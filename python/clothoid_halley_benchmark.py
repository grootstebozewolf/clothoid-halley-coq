# SPDX-FileCopyrightText: 2026 Merkator Group
# SPDX-License-Identifier: LicenseRef-Merkator-Proprietary-NoAITraining
"""
Improvement over Bertolazzi-Frego (2015) "G1 fitting with clothoids".

The canonical algorithm reduces G1 Hermite fitting to a single nonlinear
equation
    g(A) = Y_0(2A, delta-A, phi0) = 0                                  (*)
in one unknown A, solved by Newton iteration.

Idea of the improvement:
    The generalized-Fresnel routine returns all moments X_k, Y_k in one
    shot. With a tiny extra cost (nk=5 instead of 3) we get not only
    g'(A) = X_2 - X_1 but also
        g''(A) = -(Y_4 - 2 Y_3 + Y_2)
    enabling Halley's iteration (cubic convergence) for a per-step cost
    only modestly above Newton's.  Net effect: fewer iterations to reach
    machine precision.

This benchmark sweeps (phi0, phi1) over a fine grid, runs both methods
with the same initial guess (Bertolazzi-Frego's polynomial seed), and
reports iteration counts.
"""

import math
import numpy as np
from numpy.polynomial.legendre import leggauss

# 32-pt Gauss-Legendre on [0,1] -- exponential accuracy for smooth
# integrands; both methods use the same quadrature so the comparison is
# purely about iteration count.
_GL_N = 64
_xi, _wi = leggauss(_GL_N)
_T = 0.5 * (_xi + 1.0)
_W = 0.5 * _wi
_T_pow = {k: _T**k for k in range(6)}     # cache t^k


def generalized_fresnel(nk, a, b, c):
    """Return (X[0..nk-1], Y[0..nk-1]) where
       X_k = int_0^1 t^k cos((a/2) t^2 + b t + c) dt,
       Y_k = int_0^1 t^k sin((a/2) t^2 + b t + c) dt.
    """
    phi = 0.5*a*_T_pow[2] + b*_T + c
    c_ = np.cos(phi); s_ = np.sin(phi)
    X = np.empty(nk); Y = np.empty(nk)
    for k in range(nk):
        tk = _T_pow[k]
        X[k] = np.dot(_W, tk * c_)
        Y[k] = np.dot(_W, tk * s_)
    return X, Y


# Bertolazzi-Frego initial-guess polynomial (verbatim from buildClothoid.m).
_CF = (2.989696028701907, 0.716228953608281, -0.458969738821509,
       -0.502821153340377, 0.261062141752652, -0.045854475238709)

def guess_A(phi0, phi1):
    X = phi0 / math.pi
    Y = phi1 / math.pi
    xy = X*Y
    return (phi0 + phi1) * (
        _CF[0] + xy*(_CF[1] + xy*_CF[2])
        + (_CF[3] + xy*_CF[4])*(X*X + Y*Y)
        + _CF[5]*(X**4 + Y**4)
    )


# -------------- Canonical Newton (Bertolazzi-Frego) -----------------
def solve_newton(phi0, phi1, *, tol=1e-12, max_iter=100):
    delta = phi1 - phi0
    A = guess_A(phi0, phi1)
    for it in range(1, max_iter+1):
        X, Y = generalized_fresnel(3, 2*A, delta - A, phi0)
        f  = Y[0]
        fp = X[2] - X[1]
        A -= f / fp
        if abs(f) < tol:
            return A, it
    return A, max_iter


# -------------- Proposed Halley (cubic convergence) -----------------
def solve_halley(phi0, phi1, *, tol=1e-12, max_iter=100):
    delta = phi1 - phi0
    A = guess_A(phi0, phi1)
    for it in range(1, max_iter+1):
        X, Y = generalized_fresnel(5, 2*A, delta - A, phi0)
        f   = Y[0]
        fp  = X[2] - X[1]
        fpp = -(Y[4] - 2*Y[3] + Y[2])
        # Halley step: A -= 2 f f' / (2 f'^2 - f f'')
        denom = 2.0*fp*fp - f*fpp
        A -= 2.0*f*fp / denom
        if abs(f) < tol:
            return A, it
    return A, max_iter


# -------------- Verification + hypercube benchmark ------------------
def verify_pair(phi0, phi1):
    A1, it1 = solve_newton(phi0, phi1)
    A2, it2 = solve_halley(phi0, phi1)
    return A1, it1, A2, it2, abs(A1 - A2)


def hypercube_benchmark(grid=21, exclude_radius=1e-3):
    """Sweep phi0, phi1 in (-pi, pi)^2."""
    angles = np.linspace(-math.pi + 1e-4, math.pi - 1e-4, grid)
    n_newton, n_halley = [], []
    max_disc = 0.0
    n_skipped = 0
    for p0 in angles:
        for p1 in angles:
            # skip near-degenerate (both angles near 0 simultaneously - trivial straight)
            if abs(p0) < exclude_radius and abs(p1) < exclude_radius:
                n_skipped += 1
                continue
            A1, it1 = solve_newton(p0, p1)
            A2, it2 = solve_halley(p0, p1)
            n_newton.append(it1)
            n_halley.append(it2)
            max_disc = max(max_disc, abs(A1 - A2))
    return (np.array(n_newton), np.array(n_halley), max_disc, n_skipped)


if __name__ == "__main__":
    print("Sanity check on a handful of (phi0, phi1) configurations")
    print("-"*70)
    for p0, p1 in [(0.5, -0.4), (1.2, 0.8), (-1.5, 1.5), (0.01, 0.02),
                   (2.5, -2.7), (0.3, 0.3), (-0.7, 1.1)]:
        A1, it1, A2, it2, disc = verify_pair(p0, p1)
        print(f"  phi0={p0:+.3f}  phi1={p1:+.3f}  "
              f"Newton: A={A1:+.10f} ({it1} iter)   "
              f"Halley: A={A2:+.10f} ({it2} iter)   "
              f"|disc|={disc:.1e}")

    print()
    print("Hypercube benchmark on a 21x21 grid of (phi0,phi1) in (-pi,pi)^2")
    print("-"*70)
    nN, nH, max_disc, skipped = hypercube_benchmark(grid=21)
    print(f"samples: {nN.size}   skipped (both angles ~ 0): {skipped}")
    print(f"max |A_newton - A_halley| over grid: {max_disc:.2e}")
    print()
    print(f"{'method':>10s}  {'mean':>6s}  {'median':>6s}  {'p95':>6s}  "
          f"{'max':>4s}  {'>=5':>6s}  {'>=6':>6s}")
    for name, arr in [("Newton (BF)", nN), ("Halley (new)", nH)]:
        print(f"{name:>12s}  {arr.mean():>6.2f}  {np.median(arr):>6.1f}  "
              f"{np.percentile(arr,95):>6.1f}  {arr.max():>4d}  "
              f"{(arr>=5).mean()*100:>5.1f}%  {(arr>=6).mean()*100:>5.1f}%")

    # Histogram
    print("\nIteration-count histogram:")
    print(f"{'iters':>5s}  {'Newton':>8s}  {'Halley':>8s}")
    for k in range(1, max(nN.max(), nH.max()) + 1):
        print(f"{k:>5d}  {(nN == k).sum():>8d}  {(nH == k).sum():>8d}")
