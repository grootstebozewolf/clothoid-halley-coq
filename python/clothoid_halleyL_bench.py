# SPDX-FileCopyrightText: 2026 Merkator Group
# SPDX-License-Identifier: EUPL-1.2
"""
Halley vs Newton on the chord-length residual

    f(L) = L^2 (P^2 + Q^2) - d^2,    P=int cos(L psi),  Q=int sin(L psi)

with f'(L), f''(L) derived in clothoid_halleyL_derive.py (sympy-verified):

    f'(L)  = 2L (P^2+Q^2) + 2L^2 (QR - PT)
    f''(L) = 2 (P^2+Q^2) + 8L (QR - PT) + 2L^2 (R^2 + T^2 - P S2c - Q S2s)

Six moments are needed per Halley step (P, Q, R, T, S2c, S2s).  All come from
the same 32-pt Gauss-Legendre cache used for the Newton variant.

Initial guess: L_0 = d (strict lower bound, as before).
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

_GL_N = 32
_xi, _wi = leggauss(_GL_N)
_T = 0.5*(_xi + 1.0)
_W = 0.5*_wi


def moments(L, k0, k1):
    """Return P, Q, R, T, S2c, S2s on [0,1]."""
    psi = k0*_T + 0.5*(k1 - k0)*_T*_T
    c = np.cos(L*psi); s = np.sin(L*psi)
    P   = np.dot(_W, c)
    Q   = np.dot(_W, s)
    R   = np.dot(_W, psi*c)
    T   = np.dot(_W, psi*s)
    S2c = np.dot(_W, psi*psi*c)
    S2s = np.dot(_W, psi*psi*s)
    return P, Q, R, T, S2c, S2s


def solve_newton_L(P0, P1, k0, k1, *, tol=1e-13, max_iter=50):
    Cx, Cy = P1[0]-P0[0], P1[1]-P0[1]
    d2 = Cx*Cx + Cy*Cy
    d  = math.sqrt(d2)
    if d == 0: return 0.0, 0
    L = d
    for it in range(1, max_iter+1):
        P, Q, R, T, _, _ = moments(L, k0, k1)
        r2 = P*P + Q*Q
        f  = L*L*r2 - d2
        fp = 2*L*r2 + 2*L*L*(Q*R - P*T)
        if abs(f) < tol*max(d2, 1.0): return L, it
        if fp <= 0:
            L *= 1.5; continue
        L = max(L - f/fp, 0.5*L)
    return L, max_iter


def solve_halley_L(P0, P1, k0, k1, *, tol=1e-13, max_iter=50):
    Cx, Cy = P1[0]-P0[0], P1[1]-P0[1]
    d2 = Cx*Cx + Cy*Cy
    d  = math.sqrt(d2)
    if d == 0: return 0.0, 0
    L = d
    for it in range(1, max_iter+1):
        P, Q, R, T, S2c, S2s = moments(L, k0, k1)
        r2  = P*P + Q*Q
        QRPT = Q*R - P*T
        f   = L*L*r2 - d2
        fp  = 2*L*r2 + 2*L*L*QRPT
        fpp = 2*r2 + 8*L*QRPT + 2*L*L*(R*R + T*T - P*S2c - Q*S2s)
        if abs(f) < tol*max(d2, 1.0): return L, it
        denom = 2*fp*fp - f*fpp
        if abs(denom) < 1e-20 or fp <= 0:
            L *= 1.5; continue
        step = 2*f*fp / denom
        L_new = L - step
        if L_new <= 0: L_new = 0.5*L
        L = L_new
    return L, max_iter


def forward(P0, theta0, k0, k1, L):
    alpha = (k1 - k0)/L
    def h(s): return theta0 + k0*s + 0.5*alpha*s*s
    cx, _ = quad(lambda s: math.cos(h(s)), 0, L, epsabs=1e-13, epsrel=1e-13)
    cy, _ = quad(lambda s: math.sin(h(s)), 0, L, epsabs=1e-13, epsrel=1e-13)
    return (P0[0]+cx, P0[1]+cy)


def hypercube():
    """Sweep over a representative parameter cube.

    Stay on the first monotone branch of f(L) (no spiraling): keep
    k0*L, k1*L within +-pi to avoid multi-root regime.
    """
    cases = []
    for k0 in np.linspace(-1.5, 1.5, 9):
        for k1 in np.linspace(-1.5, 1.5, 9):
            for L_true in [0.5, 1.0, 2.0, 3.0]:
                if abs(k0*L_true) > math.pi or abs(k1*L_true) > math.pi:
                    continue
                # any theta0 / P0 (problem is translation- and rotation-invariant)
                P0 = (0.0, 0.0); theta0 = 0.3
                P1 = forward(P0, theta0, k0, k1, L_true)
                cases.append((P0, P1, k0, k1, L_true))
    return cases


def main():
    cases = hypercube()
    print(f"samples: {len(cases)}")
    nN, nH = [], []
    max_disc = 0.0
    worst_N = (0, 0, 0, 0)
    worst_H = (0, 0, 0, 0)
    for P0, P1, k0, k1, L_true in cases:
        LN, iN = solve_newton_L(P0, P1, k0, k1)
        LH, iH = solve_halley_L(P0, P1, k0, k1)
        nN.append(iN); nH.append(iH)
        max_disc = max(max_disc, abs(LN-LH))
        if iN > worst_N[3]: worst_N = (k0, k1, L_true, iN)
        if iH > worst_H[3]: worst_H = (k0, k1, L_true, iH)
    nN, nH = np.array(nN), np.array(nH)
    print(f"max |L_newton - L_halley| = {max_disc:.2e}")
    print(f"{'method':>14s}  mean   median   max  >=4   >=5")
    for name, arr in [("Newton", nN), ("Halley", nH)]:
        print(f"{name:>14s}  {arr.mean():.2f}   {np.median(arr):.0f}      {arr.max():2d}   "
              f"{(arr>=4).mean()*100:4.1f}%  {(arr>=5).mean()*100:4.1f}%")
    print(f"\nworst Newton: k0={worst_N[0]:.2f} k1={worst_N[1]:.2f} L={worst_N[2]} -> {worst_N[3]}")
    print(f"worst Halley: k0={worst_H[0]:.2f} k1={worst_H[1]:.2f} L={worst_H[2]} -> {worst_H[3]}")
    print()
    print("iter-count histogram:")
    print(f"{'iters':>5s}  {'Newton':>14}  {'Halley':>14}")
    for k in range(1, max(nN.max(), nH.max())+1):
        nN_k, nH_k = (nN==k).sum(), (nH==k).sum()
        pN, pH = nN_k/nN.size*100, nH_k/nH.size*100
        print(f"{k:>5d}  {nN_k:>7d} ({pN:5.1f}%)  {nH_k:>7d} ({pH:5.1f}%)")


if __name__ == "__main__":
    main()
