# SPDX-FileCopyrightText: 2026 Merkator Group
# SPDX-License-Identifier: LicenseRef-Merkator-Proprietary-NoAITraining
"""Peer review of the submitted fit_clothoid_halley solver."""
import math
import numpy as np
from scipy.special import fresnel
from scipy.integrate import quad


# ----- Code under review (verbatim from submission) ------------------------
def clothoid_pq(L, k0, k1):
    if abs(L) < 1e-14:
        return 1.0, 0.0
    a = L * k0
    b = L * (k1 - k0) / 2.0
    if abs(b) < 1e-14:
        if abs(a) < 1e-14:
            return 1.0, 0.0
        z = (np.exp(1j * a) - 1.0) / (1j * a)
        return z.real, z.imag
    u0 = a / (2.0 * b)
    phase = np.exp(1j * (-b * u0**2))
    sigma = np.sqrt(2.0 * abs(b) / np.pi)
    t1 = sigma * u0
    t2 = sigma * (u0 + 1.0)
    S2, C2 = fresnel(t2)
    S1, C1 = fresnel(t1)
    I = (C2 - C1) + 1j * (S2 - S1)
    I /= sigma
    if b < 0.0:
        I = np.conj(I)
    I *= phase
    return I.real, I.imag


def fit_clothoid_halley(P0, P1, k0, k1, L_target, eps=1e-6, tol=1e-14, maxit=4):
    dx = P1[0] - P0[0]
    dy = P1[1] - P0[1]
    d2 = dx*dx + dy*dy
    d = math.hypot(dx, dy)
    if d < tol:
        return 0.0, 0.0
    L = L_target
    Lmin = max(0.0, L_target - eps)
    Lmax = L_target + eps
    iters_used = 0
    for _ in range(maxit):
        iters_used += 1
        P, Q = clothoid_pq(L, k0, k1)
        r2 = P*P + Q*Q
        f = L*L * r2 - d2
        if abs(f) < tol * max(d2, 1.0):
            break
        h = 1e-9
        Pp, Qp = clothoid_pq(L + h, k0, k1)
        r2p = Pp*Pp + Qp*Qp
        fp = ((L + h)**2 * r2p - L**2 * r2) / h
        if abs(fp) < 1e-20:
            L = max(Lmin, min(Lmax, L * d / (L * math.sqrt(r2) + 1e-20)))
            continue
        step = f / fp                       # <<< this is Newton, not Halley
        L -= step
        L = max(Lmin, min(Lmax, L))
    else:
        return None, None, iters_used
    P, Q = clothoid_pq(L, k0, k1)
    theta0 = math.atan2(dy, dx) - math.atan2(Q, P)
    return L, theta0, iters_used


# ----- Reference forward integrator (ground truth) -------------------------
def forward(P0, theta0, k0, k1, L):
    alpha = (k1 - k0) / L
    def h(s): return theta0 + k0*s + 0.5*alpha*s*s
    cx, _ = quad(lambda s: math.cos(h(s)), 0, L, epsabs=1e-14, epsrel=1e-14)
    cy, _ = quad(lambda s: math.sin(h(s)), 0, L, epsabs=1e-14, epsrel=1e-14)
    return (P0[0]+cx, P0[1]+cy)


# ----- A. Verify clothoid_pq against direct quadrature ---------------------
print("="*72)
print("A. clothoid_pq vs direct quadrature (sanity check on the Fresnel core)")
print("="*72)
for L, k0, k1 in [(4.0, 0.0, 0.5), (3.0, 0.7, -0.7), (2.5, 0.4, 0.4),
                  (7.0, 0.0, 0.0), (1.5, 2.0, -1.0), (50.0, 0.05, 0.1)]:
    P, Q = clothoid_pq(L, k0, k1)
    Pref, _ = quad(lambda t: math.cos(L*(k0*t + 0.5*(k1-k0)*t*t)), 0, 1,
                   epsabs=1e-13, epsrel=1e-13)
    Qref, _ = quad(lambda t: math.sin(L*(k0*t + 0.5*(k1-k0)*t*t)), 0, 1,
                   epsabs=1e-13, epsrel=1e-13)
    print(f"  L={L:>4}  k0={k0:+.2f} k1={k1:+.2f}  "
          f"|dP|={abs(P-Pref):.1e}  |dQ|={abs(Q-Qref):.1e}")


# ----- B. CLAIM #1: "Halley" -> is it actually Halley? ---------------------
print()
print("="*72)
print("B. The step formula is f/fp, i.e. Newton.  Halley would be 2*f*fp/(2*fp^2 - f*fpp).")
print("   No fpp is ever computed.  The function name is misleading.")
print("="*72)


# ----- C. CLAIM #2: "machine-precision recovery on the original cases" -----
print()
print("="*72)
print("C. Test claim: passes the original test suite at machine precision")
print("="*72)
# The original suite (clothoid_newton_verify.py) gave (P0, P1, k0, k1) ONLY -
# no L_target.  So the realistic comparison is fit_clothoid_halley with
# L_target equal to... what?  The submission doesn't say.  The claim
# "Recovered L*=4.0 ... recovers exactly" only works if L_target is already
# very close to L_true.  Test both regimes.

cases = [
    ("k0=0   k1=0.5  L=4",   (0,0), 0.3,  0.0,  0.5, 4.0),
    ("S-curve k0=0.7 k1=-0.7 L=3",  (0,0), 1.0,  0.7, -0.7, 3.0),
    ("pure arc k0=k1=0.4 L=2.5", (0,0), 0.0,  0.4,  0.4, 2.5),
    ("strong curve k0=2 k1=-1 L=1.5",(3,-4),-2.0, 2.0, -1.0, 1.5),
    ("long shallow k0=0.05 k1=0.1 L=50",(0,0),0.0, 0.05, 0.1, 50.0),
]

print("\n  C1.  L_target = L_true (the easy/dishonest test)")
print("       this should always succeed in 1-2 iters")
print(f"  {'case':<36s}  {'L_rec':>8s}  {'|dL|':>9s}  iters")
for name, P0, t0, k0, k1, L_true in cases:
    P1 = forward(P0, t0, k0, k1, L_true)
    res = fit_clothoid_halley(P0, P1, k0, k1, L_target=L_true, eps=1e-6, maxit=4)
    L, th, it = res if res[0] is not None else (None, None, "fail")
    if L is None:
        print(f"  {name:<36s}  {'fail':>8s}")
    else:
        print(f"  {name:<36s}  {L:>8.4f}  {abs(L-L_true):>9.1e}  {it}")

print("\n  C2.  L_target = chord length d (realistic blind seed used in the")
print("       original solver -- L_target should NOT be the answer)")
print(f"  {'case':<36s}  {'L_rec':>10s}  {'|dL|':>9s}  iters  band")
for name, P0, t0, k0, k1, L_true in cases:
    P1 = forward(P0, t0, k0, k1, L_true)
    d = math.hypot(P1[0]-P0[0], P1[1]-P0[1])
    # narrow band first
    res = fit_clothoid_halley(P0, P1, k0, k1, L_target=d, eps=1e-6, maxit=4)
    L, th, it = (res[0], res[1], res[2]) if res[0] is not None else (None, None, "fail")
    label = "narrow"
    print(f"  {name:<36s}  {str(L):>10s}  "
          f"{'' if L is None else f'{abs(L-L_true):.1e}':>9s}  {it}  {label}")
    # wide band
    res = fit_clothoid_halley(P0, P1, k0, k1, L_target=d, eps=L_true*2, maxit=4)
    L, th, it = (res[0], res[1], res[2]) if res[0] is not None else (None, None, "fail")
    label = "eps=2L"
    print(f"  {'':<36s}  {str(L):>10s}  "
          f"{'' if L is None else f'{abs(L-L_true):.1e}':>9s}  {it}  {label}")


# ----- D. CLAIM #3: "FD with h=1e-9 is analytic-quality" -------------------
print()
print("="*72)
print("D. FD derivative quality vs analytic (h=1e-9)")
print("="*72)
def analytic_fp(L, k0, k1, d2):
    # exact: fp = 2L(P^2+Q^2) + 2L^2(Q*R - P*T), with R,T = first moments
    from numpy.polynomial.legendre import leggauss
    xi, wi = leggauss(64)
    T_ = 0.5*(xi + 1.0); W = 0.5*wi
    psi = k0*T_ + 0.5*(k1-k0)*T_*T_
    c = np.cos(L*psi); s = np.sin(L*psi)
    P = np.dot(W, c); Q = np.dot(W, s)
    R = np.dot(W, psi*c); T = np.dot(W, psi*s)
    return 2*L*(P*P + Q*Q) + 2*L*L*(Q*R - P*T)

def fd_fp(L, k0, k1, d2, h):
    P,  Q  = clothoid_pq(L,   k0, k1)
    Pp, Qp = clothoid_pq(L+h, k0, k1)
    return ((L+h)**2 * (Pp*Pp+Qp*Qp) - L*L * (P*P+Q*Q)) / h

print(f"  {'(L, k0, k1)':<24s}  {'analytic':>14}  {'FD h=1e-9':>14}  "
      f"{'FD h=1e-6':>14}  {'rel err 1e-9':>12}")
for L, k0, k1 in [(4.0, 0.0, 0.5), (50.0, 0.05, 0.1), (1.5, 2.0, -1.0)]:
    a = analytic_fp(L, k0, k1, 0)
    f1 = fd_fp(L, k0, k1, 0, 1e-9)
    f2 = fd_fp(L, k0, k1, 0, 1e-6)
    print(f"  L={L:<5} k0={k0:<5} k1={k1:<5}  {a:>14.6e}  {f1:>14.6e}  "
          f"{f2:>14.6e}  {abs(a-f1)/abs(a):>12.1e}")


# ----- E. CLAIM #4: "long shallow ... Recovered L*=50, ΔL=0" ---------------
print()
print("="*72)
print("E. The 'long shallow' case has multiple roots of f(L).  Re-examine.")
print("="*72)
# Find any L >= d where f(L)=0 by sweeping
P0=(0,0); k0=0.05; k1=0.1
P1 = forward(P0, 0.0, k0, k1, 50.0)
dx, dy = P1[0]-P0[0], P1[1]-P0[1]
d2 = dx*dx + dy*dy
def f_of_L(L):
    P, Q = clothoid_pq(L, k0, k1)
    return L*L*(P*P+Q*Q) - d2

Ls = np.linspace(1.0, 60.0, 600)
sign_changes = []
prev = f_of_L(Ls[0])
for L in Ls[1:]:
    cur = f_of_L(L)
    if prev*cur < 0:
        sign_changes.append(L)
    prev = cur
print(f"  sign changes of f(L) on [1,60]: {sign_changes}")
print(f"  chord d = sqrt(d2) = {math.sqrt(d2):.4f}")
print(f"  -> f has multiple roots; only one matches L_true=50.")
print(f"  with L_target=50 (= L_true), Newton terminates immediately,")
print(f"  which is what produces the |dL|=0 row in the submission's table.")
print(f"  This isn't a test of the solver -- it's a test of f(50) ~ 0.")
