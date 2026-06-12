# SPDX-FileCopyrightText: 2026 Merkator Group
# SPDX-License-Identifier: EUPL-1.2
"""
Symbolic derivation of f, f', f'' for the chord-length residual

    f(L) = L^2 (P^2 + Q^2) - d^2

with P, Q defined as functions of L via

    psi(tau)        = k0*tau + (k1 - k0)*tau^2 / 2     (independent of L)
    P(L) = int_0^1 cos(L*psi(tau)) dtau
    Q(L) = int_0^1 sin(L*psi(tau)) dtau

This is needed for Halley on L (the second problem you asked about, where
P0, P1, k0, k1 are fixed and L is the unknown).

Goal: produce f''(L) in a form expressible in standard moments
    R   = int psi   cos(L psi),   T   = int psi   sin(L psi)
    S2c = int psi^2 cos(L psi),   S2s = int psi^2 sin(L psi)
and compare against the user-supplied formula
    f''(L) ?=  2 r^2 + 8 L (Q R - P T) + 2 L^2 ( Q*S2c - P*S2s - (R^2 + T^2) )
"""
import sympy as sp

L, tau, k0, k1, d = sp.symbols('L tau kappa_0 kappa_1 d', real=True, positive=True)
psi = k0*tau + (k1 - k0)*tau**2/2

# Build P, Q, R, T, S2c, S2s as integrals in sympy
P   = sp.Function('P')(L)
Q   = sp.Function('Q')(L)
# We will compute their L-derivatives by differentiating the integrand under
# the (constant) integration limits and re-expressing in moments.

# Differentiate integrand for each moment
def D(integrand): return sp.diff(integrand, L)

cosLpsi = sp.cos(L*psi)
sinLpsi = sp.sin(L*psi)

# Direct moment definitions
P_int   = sp.integrate(cosLpsi,      (tau, 0, 1))   # = P(L)
Q_int   = sp.integrate(sinLpsi,      (tau, 0, 1))   # = Q(L)

# Derivatives expressed as moments
# dP/dL    = int -psi sin(L psi) = -T
# dQ/dL    = int  psi cos(L psi) =  R
# dR/dL    = int -psi^2 sin(L psi) = -S2s
# dT/dL    = int  psi^2 cos(L psi) =  S2c
# Build f, f', f'' fully symbolically (no integrals -- name placeholders)
Ps, Qs, Rs, Ts, S2cs, S2ss = sp.symbols('P Q R T S2c S2s', real=True)
# Encode the rules d/dL of each variable as a sympy substitution-aware chain
# We will manually take f''(L) using product rule, then substitute the
# derivative identities above.

f      = L**2 * (Ps**2 + Qs**2) - d**2
# Build f'(L) by hand using P'=-T, Q'=R
fp     = 2*L*(Ps**2 + Qs**2) + 2*L**2*(Qs*Rs - Ps*Ts)
print("f'(L)  =", sp.simplify(fp))

# d/dL of fp:
#   d/dL [2L(P^2+Q^2)]      = 2(P^2+Q^2) + 2L*(2P*(-T)+2Q*R) = 2 r^2 + 4L(QR-PT)
#   d/dL [2L^2(QR - PT)]    = 4L(QR - PT) + 2L^2 d/dL[QR - PT]
#   d/dL[QR-PT]             = (R)*R + Q*(-S2s) - (-T)*T - P*S2c
#                           = R^2 + T^2 - Q*S2s - P*S2c
# So
#   f''(L) = 2 r^2 + 8L(QR - PT) + 2L^2 ( R^2 + T^2 - P*S2c - Q*S2s )

fpp_expected = (2*(Ps**2 + Qs**2)
                + 8*L*(Qs*Rs - Ps*Ts)
                + 2*L**2*(Rs**2 + Ts**2 - Ps*S2cs - Qs*S2ss))
print("f''(L) =", fpp_expected)

# User-supplied formula:
fpp_user = (2*(Ps**2 + Qs**2)
            + 8*L*(Qs*Rs - Ps*Ts)
            + 2*L**2*(Qs*S2cs - Ps*S2ss - (Rs**2 + Ts**2)))
print("\nUser formula f''(L) =", fpp_user)
diff = sp.expand(fpp_expected - fpp_user)
print("\n(correct) - (user) =", diff)
print("Same formula?", diff == 0)

# Independent cross-check via direct symbolic differentiation of f(L) where
# P and Q are real integrals.  Use a concrete (k0, k1) so sympy can simplify.
print("\n--- Numerical cross-check at a concrete (k0,k1,L) ---")
import sympy as sp
import math
import scipy.integrate as si

def num_moments(L_val, k0_val, k1_val):
    def psi(t): return k0_val*t + (k1_val-k0_val)*t*t/2
    def integ(g): return si.quad(g, 0, 1, epsabs=1e-13, epsrel=1e-13)[0]
    P_v   = integ(lambda t: math.cos(L_val*psi(t)))
    Q_v   = integ(lambda t: math.sin(L_val*psi(t)))
    R_v   = integ(lambda t: psi(t)*math.cos(L_val*psi(t)))
    T_v   = integ(lambda t: psi(t)*math.sin(L_val*psi(t)))
    S2c_v = integ(lambda t: psi(t)**2*math.cos(L_val*psi(t)))
    S2s_v = integ(lambda t: psi(t)**2*math.sin(L_val*psi(t)))
    return P_v, Q_v, R_v, T_v, S2c_v, S2s_v

def num_f(L_val, k0_val, k1_val, d_val):
    P_v, Q_v, *_ = num_moments(L_val, k0_val, k1_val)
    return L_val**2*(P_v**2 + Q_v**2) - d_val**2

# Build f, f', f'' three ways and compare:
# (1) finite differences of f
# (2) the candidate expected formula evaluated on moments
# (3) the user formula evaluated on moments
import numpy as np
for k0v, k1v, Lv in [(0.3, -0.4, 2.5), (1.0, 0.5, 3.0), (-0.7, 0.7, 4.0)]:
    P_v, Q_v, R_v, T_v, S2c_v, S2s_v = num_moments(Lv, k0v, k1v)
    d_v = math.sqrt(Lv**2*(P_v**2 + Q_v**2))   # so f(L*) = 0; that L* = Lv

    # finite-diff f''
    h = 1e-4
    fpp_fd = (num_f(Lv+h, k0v, k1v, d_v)
              - 2*num_f(Lv, k0v, k1v, d_v)
              + num_f(Lv-h, k0v, k1v, d_v)) / h**2

    fpp_mine = (2*(P_v**2 + Q_v**2)
                + 8*Lv*(Q_v*R_v - P_v*T_v)
                + 2*Lv**2*(R_v**2 + T_v**2 - P_v*S2c_v - Q_v*S2s_v))
    fpp_user_ = (2*(P_v**2 + Q_v**2)
                 + 8*Lv*(Q_v*R_v - P_v*T_v)
                 + 2*Lv**2*(Q_v*S2c_v - P_v*S2s_v - (R_v**2 + T_v**2)))
    print(f"k0={k0v:+.2f} k1={k1v:+.2f} L={Lv}:")
    print(f"  f''(fd)    = {fpp_fd:+.10f}")
    print(f"  f''(mine)  = {fpp_mine:+.10f}   diff_fd = {abs(fpp_mine-fpp_fd):.1e}")
    print(f"  f''(user)  = {fpp_user_:+.10f}   diff_fd = {abs(fpp_user_-fpp_fd):.1e}")
