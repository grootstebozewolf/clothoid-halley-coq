"""
Deterministic clothoid algorithm derivation using sympy.

Inputs:  P0 = (x0,y0), P1 = (x1,y1), kappa0, kappa1, L
Output:  parametric curve P(s), 0<=s<=L, single-pass (no Newton iteration).

Canonical algorithms for the G1 Hermite clothoid (given P0, theta0, P1, theta1)
solve a transcendental equation by Newton iteration. Here we change the input set:
the heading theta0 is unknown but is recoverable in closed form from the chord
when (kappa0, kappa1, L) are given. The remaining integrals are pure Fresnel
integrals whose argument depends only on (kappa0, kappa1, L), so no nonlinear
solver is needed.
"""

import sympy as sp

s, t, L, k0, k1, theta0 = sp.symbols('s t L kappa_0 kappa_1 theta_0', real=True)
x0, y0, x1, y1 = sp.symbols('x_0 y_0 x_1 y_1', real=True)

# --- 1. Clothoid definition: curvature linear in arclength -------------------
alpha = (k1 - k0) / L                                  # curvature rate
kappa = k0 + alpha * s
theta = theta0 + sp.integrate(k0 + alpha*t, (t, 0, s)) # heading
theta = sp.simplify(theta)
print("kappa(s) =", kappa)
print("theta(s) =", theta)

# --- 2. Position integrals ---------------------------------------------------
#   x(s) = x0 + Integral cos(theta(t)) dt
#   y(s) = y0 + Integral sin(theta(t)) dt
# Factor out theta0 so the integrand depends only on (k0,k1,L):
phi = theta - theta0
print("phi(s) = theta(s) - theta0 =", sp.simplify(phi))

# Define the "base" Fresnel-like integrals (no theta0):
Ic = sp.Function('I_c')(s)     # = integral_0^s cos(phi(t)) dt
Is = sp.Function('I_s')(s)     # = integral_0^s sin(phi(t)) dt

# Then by angle-addition:
#   x(s) - x0 = cos(theta0)*Ic(s) - sin(theta0)*Is(s)
#   y(s) - y0 = sin(theta0)*Ic(s) + cos(theta0)*Is(s)

# --- 3. Reduce phi(t) = k0*t + alpha*t^2/2 to a pure Fresnel argument --------
# Complete the square: phi(t) = (alpha/2)*(t + k0/alpha)^2 - k0^2/(2*alpha)
# Substitute u = sqrt(|alpha|/pi) * (t + k0/alpha)  =>  pi*u^2/2 = (alpha/2)*(t+k0/alpha)^2
u = sp.symbols('u', real=True)
t_of_u = sp.sqrt(sp.pi/sp.Abs(alpha)) * u - k0/alpha
# show the shift constant
shift = -k0**2/(2*alpha)
print("shift constant gamma = -k0^2/(2*alpha) =", sp.simplify(shift))
# After substitution, phi(t) = pi*u^2/2 * sign(alpha) + shift, dt = sqrt(pi/|alpha|) du

# Limits: when t=0, u0 = k0/sqrt(pi*|alpha|)*sign(alpha)... handle sign carefully.
# For the symbolic derivation assume alpha>0 (the alpha<0 case mirrors via sign).
u_lo = sp.sqrt(sp.Abs(alpha)/sp.pi) * (0 + k0/alpha)
u_hi = sp.sqrt(sp.Abs(alpha)/sp.pi) * (s + k0/alpha)
print("u_lo =", sp.simplify(u_lo))
print("u_hi =", sp.simplify(u_hi))

# Standard Fresnel integrals:
#   C(u) = integral_0^u cos(pi*v^2/2) dv,   S(u) = integral_0^u sin(pi*v^2/2) dv
C = sp.Function('C')   # Fresnel C
S = sp.Function('S')   # Fresnel S

# So (assuming alpha>0):
#   Ic(s) = sqrt(pi/alpha) * [ cos(shift)*(C(u_hi)-C(u_lo)) - sin(shift)*(S(u_hi)-S(u_lo)) ]
#   Is(s) = sqrt(pi/alpha) * [ cos(shift)*(S(u_hi)-S(u_lo)) + sin(shift)*(C(u_hi)-C(u_lo)) ]

# --- 4. Closed-form recovery of theta0 from the chord ------------------------
# With Cx = x1-x0, Cy = y1-y0 and the *terminal* integrals A=Ic(L), B=Is(L):
A, B, Cx, Cy = sp.symbols('A B C_x C_y', real=True)
sol_theta0 = sp.atan2(Cy, Cx) - sp.atan2(B, A)
print("theta0 =", sol_theta0)

# Consistency check: |chord|^2 must equal A^2 + B^2.
print("consistency: Cx^2 + Cy^2 ==", sp.simplify(A**2 + B**2), "(must hold)")

# --- 5. Pulling it all together: deterministic point evaluator ---------------
print("\n=== Deterministic clothoid point P(s) ===")
P_x = x0 + sp.cos(theta0)*Ic - sp.sin(theta0)*Is
P_y = y0 + sp.sin(theta0)*Ic + sp.cos(theta0)*Is
print("x(s) = ", P_x)
print("y(s) = ", P_y)
