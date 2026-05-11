"""
Newton's method on the length L for the clothoid connecting P0->P1 with
prescribed (kappa_0, kappa_1).

Iteration:
    L_{n+1} = L_n - f(L_n)/f'(L_n)
    f(L)   = L^2 (P^2 + Q^2) - d^2          (squared chord-length residual)
    f'(L)  = 2L (P^2 + Q^2) + 2 L^2 (Q*R - P*T)
    where, with psi(tau) = k0*tau + (k1-k0)*tau^2/2,
        P = int_0^1 cos(L*psi(tau)) dtau
        Q = int_0^1 sin(L*psi(tau)) dtau
        R = int_0^1 psi(tau)*cos(L*psi(tau)) dtau
        T = int_0^1 psi(tau)*sin(L*psi(tau)) dtau

P, Q reduce to standard Fresnel integrals (see earlier derivation).
R, T are first moments; here we use fixed Gauss-Legendre (32-pt), which gives
constant per-iteration cost. Quadratic Newton convergence => ~log2(N) iterations
to reach N bits of precision.
"""

import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

# 32-point Gauss-Legendre on [0,1] — fixed cost per Newton step
_GL_N = 32
_xi, _wi = leggauss(_GL_N)
_nodes  = 0.5*(_xi + 1.0)
_weights = 0.5*_wi


def _PQRT(L, k0, k1):
    psi = k0*_nodes + 0.5*(k1 - k0)*_nodes**2
    c = np.cos(L*psi); s = np.sin(L*psi)
    P = np.dot(_weights, c)
    Q = np.dot(_weights, s)
    R = np.dot(_weights, psi*c)
    T = np.dot(_weights, psi*s)
    return P, Q, R, T


def solve_length(P0, P1, k0, k1, *, tol=1e-13, max_iter=50, verbose=False):
    Cx, Cy = P1[0]-P0[0], P1[1]-P0[1]
    d2 = Cx*Cx + Cy*Cy
    d  = math.sqrt(d2)
    if d == 0:
        return 0.0, 0, []
    L = d                                   # strict lower bound -- safe init
    history = []
    for it in range(1, max_iter+1):
        P, Q, R, T = _PQRT(L, k0, k1)
        f  = L*L*(P*P + Q*Q) - d2
        fp = 2*L*(P*P + Q*Q) + 2*L*L*(Q*R - P*T)
        history.append((it, L, f))
        if abs(f) < tol*max(d2, 1.0):
            break
        if fp <= 0:
            # Defensive: f starts at f(d) <= 0 and we want to grow L.
            # Take a damped expansion step instead of dividing by a bad slope.
            L *= 1.5
            continue
        step = f/fp
        # Newton step; guard against overshooting past d (f would become +infty
        # well before that for non-spiraling curves).
        L_new = L - step
        if L_new <= 0:
            L_new = 0.5*L
        L = L_new
        if verbose:
            print(f"  it={it:2d}  L={L:.15f}  f={f:+.3e}")
    return L, it, history


# ---------- Reference: build the curve forward and use as ground truth ------
def forward_clothoid(P0, theta0, k0, k1, L):
    alpha = (k1 - k0)/L
    def heading(s): return theta0 + k0*s + 0.5*alpha*s*s
    cx, _ = quad(lambda s: math.cos(heading(s)), 0, L, epsabs=1e-14, epsrel=1e-14)
    cy, _ = quad(lambda s: math.sin(heading(s)), 0, L, epsabs=1e-14, epsrel=1e-14)
    return (P0[0]+cx, P0[1]+cy)


def run_case(name, P0, theta0_true, k0, k1, L_true):
    P1 = forward_clothoid(P0, theta0_true, k0, k1, L_true)
    L_rec, iters, hist = solve_length(P0, P1, k0, k1, verbose=False)
    rel = abs(L_rec - L_true)/L_true
    print(f"{name:38s}  L*={L_true:.6f}  L_hat={L_rec:.15f}  "
          f"|err|/L = {rel:.2e}  iters = {iters}")
    # Show how the residual collapses (quadratic Newton => digit-doubling)
    for it, L, f in hist[:8]:
        print(f"    it {it:2d}: L={L:.15f}  f={f:+.3e}")
    print()


print("Newton on L: P0, P1, k0, k1 given; L is the unknown")
print("="*78)
run_case("clothoid k0=0   k1=0.5  L=4",  (0,0), 0.3,  0.0,  0.5, 4.0)
run_case("clothoid k0=-0.3 k1=0.7 L=5", (1,2), -0.5,-0.3,  0.7, 5.0)
run_case("S-curve k0=0.7 k1=-0.7 L=3",  (0,0), 1.0,  0.7, -0.7, 3.0)
run_case("pure arc k0=k1=0.4 L=2.5",    (0,0), 0.0,  0.4,  0.4, 2.5)
run_case("straight k0=k1=0 L=7",        (0,0), 1.2,  0.0,  0.0, 7.0)
run_case("strong curve k0=2 k1=-1 L=1.5",(3,-4),-2.0, 2.0, -1.0, 1.5)
run_case("long shallow k0=0.05 k1=0.1 L=50",(0,0),0.0, 0.05, 0.1, 50.0)
