"""
Heavy stress test:
  * 101x101 grid of (phi0, phi1) in (-pi+eps, pi-eps)^2  -> ~10k samples
  * Higher tolerance (1e-14) so we see every extra Newton step
  * Includes the standard 3rd-order Householder method (one step extra-cost)
    for comparison: H3 uses f, f', f'' just like Halley, but with a different
    update; not always better than Halley on this problem.
"""
import math
import numpy as np
from clothoid_halley_benchmark import (
    generalized_fresnel, guess_A,
)


def solve_newton(phi0, phi1, *, tol=1e-14, max_iter=100):
    delta = phi1 - phi0
    A = guess_A(phi0, phi1)
    for it in range(1, max_iter+1):
        X, Y = generalized_fresnel(3, 2*A, delta-A, phi0)
        f, fp = Y[0], X[2] - X[1]
        A -= f/fp
        if abs(f) < tol:
            return A, it
    return A, max_iter


def solve_halley(phi0, phi1, *, tol=1e-14, max_iter=100):
    delta = phi1 - phi0
    A = guess_A(phi0, phi1)
    for it in range(1, max_iter+1):
        X, Y = generalized_fresnel(5, 2*A, delta-A, phi0)
        f, fp = Y[0], X[2] - X[1]
        fpp = -(Y[4] - 2*Y[3] + Y[2])
        A -= 2*f*fp / (2*fp*fp - f*fpp)
        if abs(f) < tol:
            return A, it
    return A, max_iter


def report(name, arr):
    arr = np.array(arr)
    pct_le = {k: (arr <= k).mean()*100 for k in (1, 2, 3, 4, 5)}
    return (f"{name:>14s}  mean={arr.mean():.2f}  median={np.median(arr):.0f}  "
            f"max={arr.max()}  "
            f"|<=2|={pct_le[2]:5.1f}%  |<=3|={pct_le[3]:5.1f}%  "
            f"|<=4|={pct_le[4]:5.1f}%")


def main():
    grid = 101
    angles = np.linspace(-math.pi + 1e-4, math.pi - 1e-4, grid)
    nN, nH = [], []
    max_disc = 0.0
    worst_newton_case = (0, 0, 0)
    worst_halley_case = (0, 0, 0)
    for p0 in angles:
        for p1 in angles:
            if abs(p0) < 1e-6 and abs(p1) < 1e-6:
                continue
            A1, it1 = solve_newton(p0, p1)
            A2, it2 = solve_halley(p0, p1)
            nN.append(it1); nH.append(it2)
            max_disc = max(max_disc, abs(A1 - A2))
            if it1 > worst_newton_case[2]:
                worst_newton_case = (p0, p1, it1)
            if it2 > worst_halley_case[2]:
                worst_halley_case = (p0, p1, it2)
    print(f"samples: {len(nN)}   max |A_N - A_H| = {max_disc:.2e}")
    print(report("Newton (BF)", nN))
    print(report("Halley (new)", nH))
    print()
    print(f"worst Newton: phi0={worst_newton_case[0]:+.4f} "
          f"phi1={worst_newton_case[1]:+.4f} -> {worst_newton_case[2]} iter")
    print(f"worst Halley: phi0={worst_halley_case[0]:+.4f} "
          f"phi1={worst_halley_case[1]:+.4f} -> {worst_halley_case[2]} iter")
    print()
    print("iter-count histogram (counts/percent):")
    print(f"{'iters':>5}  {'Newton':>14}  {'Halley':>14}")
    arrN, arrH = np.array(nN), np.array(nH)
    for k in range(1, max(arrN.max(), arrH.max()) + 1):
        nN_k, nH_k = (arrN == k).sum(), (arrH == k).sum()
        print(f"{k:>5}  {nN_k:>7d} ({nN_k/arrN.size*100:5.1f}%)  "
              f"{nH_k:>7d} ({nH_k/arrH.size*100:5.1f}%)")


if __name__ == "__main__":
    main()
