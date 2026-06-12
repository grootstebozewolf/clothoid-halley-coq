# SPDX-FileCopyrightText: 2026 Merkator Group
# SPDX-License-Identifier: EUPL-1.2
"""Derive f(L) and f'(L) for Newton iteration on the unknown length L."""
import sympy as sp

L, tau, k0, k1, d = sp.symbols('L tau kappa_0 kappa_1 d', real=True, positive=True)
psi = k0*tau + (k1 - k0)*tau**2/2

# Base integrals reparametrised by tau = s/L
P = sp.Integral(sp.cos(L*psi), (tau, 0, 1))
Q = sp.Integral(sp.sin(L*psi), (tau, 0, 1))
A = L*P
B = L*Q
f = A**2 + B**2 - d**2

# Differentiate by Leibniz under the integral. P and Q have no L-dependent
# limits, so we just differentiate the integrand.
dP = sp.diff(sp.cos(L*psi), L)              # = -psi*sin(L*psi)
dQ = sp.diff(sp.sin(L*psi), L)              # =  psi*cos(L*psi)

# Moment integrals (cost = same as Fresnel; see notes below)
T = sp.Integral(psi*sp.sin(L*psi), (tau, 0, 1))   # P'(L) = -T
R = sp.Integral(psi*sp.cos(L*psi), (tau, 0, 1))   # Q'(L) =  R

# d/dL [L*P] = P + L*(-T), so
dA = P - L*T
dB = Q + L*R
df = 2*A*dA + 2*B*dB
df = sp.expand(df)
print("f'(L) = 2*A*(P - L*T) + 2*B*(Q + L*R)")
print("      = 2L*(P^2 + Q^2) + 2*L^2*(Q*R - P*T)")

# Verify the algebraic simplification
df_simplified = 2*L*(P**2 + Q**2) + 2*L**2*(Q*R - P*T)
delta = sp.simplify(df - df_simplified)
print("symbolic check (should be 0):", delta)

# Closed form for the moments R, T in terms of Fresnel:
#   psi(tau) = (kappa_1 - kappa_0)/2 * (tau + kappa_0/(kappa_1-kappa_0))^2  - kappa_0^2/(2(kappa_1-kappa_0))
# so psi*cos(L*psi) and psi*sin(L*psi) become first-moment Fresnel integrals,
# which reduce to a Fresnel pair plus a boundary term via integration by parts:
#   d/dtau[sin(L*psi)] = L*psi'(tau)*cos(L*psi)
# so psi'(tau)*cos(L*psi) = (1/L)*d/dtau[sin(L*psi)]
# Note psi'(tau) = k0 + (k1-k0)*tau, while we have psi(tau) in the moment.
# Write psi = (1/(k1-k0))*(psi'*(k1-k0)*tau/... ) -- easier: linear combo
#   psi(tau)            = a*psi'(tau) + b*tau + c    with suitable a,b,c?
# Actually psi(tau) = tau*(k0 + (k1-k0)*tau/2) = tau*psi'(tau) - (k1-k0)*tau^2/2
# Hmm. Practical implementation just evaluates R, T by direct quadrature on
# [0,1]; per-iteration cost is constant, so this does not affect the O(log N)
# Newton rate.

# Initial guess analysis: f(d) = A(d)^2 + B(d)^2 - d^2.
# For a clothoid that does not loop, |A^2+B^2| <= L^2 with equality only for a
# straight line, so f(d) <= 0 with equality iff k0=k1=0. So L=d is always a
# strict lower bound and a safe starting point for Newton.
print("Initial guess: L_0 = d (chord length) is a strict lower bound.")
