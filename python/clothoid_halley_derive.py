# SPDX-FileCopyrightText: 2026 Merkator Group
# SPDX-License-Identifier: LicenseRef-Merkator-Proprietary-NoAITraining
"""Derive f, f', f'' symbolically for the Bertolazzi-Frego g(A) equation."""
import sympy as sp

A, t, delta, phi0 = sp.symbols('A t delta phi_0', real=True)
phi = A*t**2 + (delta - A)*t + phi0
f_integrand = sp.sin(phi)

f   = sp.Integral(f_integrand, (t, 0, 1))
fA  = sp.diff(f_integrand, A)
fAA = sp.diff(f_integrand, A, 2)

print("f(A)   integrand = sin(A t^2 + (delta-A) t + phi0)")
print("f'(A)  integrand =", sp.simplify(fA))
print("f''(A) integrand =", sp.simplify(fAA))

# Reduce f'(A) and f''(A) to moments X_k, Y_k of the generalized Fresnel:
#    X_k = int_0^1 t^k cos(phi) dt,   Y_k = int_0^1 t^k sin(phi) dt
# By inspection:
#   f'   = (t^2 - t)*cos(phi)         => X_2 - X_1
#   f''  = -(t^2 - t)^2 * sin(phi)    => -(Y_4 - 2 Y_3 + Y_2)
print("\nReduction:")
print("  f'(A)  = X_2 - X_1                 (matches Bertolazzi-Frego)")
print("  f''(A) = -(Y_4 - 2 Y_3 + Y_2)      (new -- needed for Halley)")

# Verify the f'' expansion of (t^2-t)^2
e = sp.expand((t**2 - t)**2)
print("  (t^2 - t)^2 =", e, "  => moments Y_4 - 2 Y_3 + Y_2")
