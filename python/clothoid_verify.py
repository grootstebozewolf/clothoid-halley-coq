"""
Numerical verification of the deterministic clothoid algorithm.

Algorithm (input: P0, P1, k0, k1, L):

  1. alpha = (k1 - k0) / L                              # curvature rate
  2. Compute base integrals on [0, L] (no theta0):
        A = Ic(L) = ∫_0^L cos(k0*t + alpha*t^2/2) dt
        B = Is(L) = ∫_0^L sin(k0*t + alpha*t^2/2) dt
     For |alpha| > 0 this is a closed-form Fresnel expression:
        gamma = -k0^2 / (2*alpha)
        sigma = sign(alpha)
        scale = sqrt(pi / |alpha|)
        u_lo  = sqrt(|alpha|/pi) * (k0/alpha)
        u_hi  = sqrt(|alpha|/pi) * (L + k0/alpha)
        dC    = C(u_hi) - C(u_lo)
        dS    = sigma*(S(u_hi) - S(u_lo))
        A     = scale * ( cos(gamma)*dC - sin(gamma)*dS )
        B     = scale * ( sin(gamma)*dC + cos(gamma)*dS )
     For alpha == 0:
        if k0 == 0:  A = L,                     B = 0                 (straight)
        else:        A = sin(k0*L)/k0,          B = (1 - cos(k0*L))/k0 (arc)
  3. theta0 = atan2(Cy, Cx) - atan2(B, A)          # closed form
  4. For any s in [0, L], compute (Ic(s), Is(s)) with the same formula but with
     u_hi using s instead of L, then
        x(s) = x0 + cos(theta0)*Ic(s) - sin(theta0)*Is(s)
        y(s) = y0 + sin(theta0)*Ic(s) + cos(theta0)*Is(s)

This is deterministic: no Newton iteration. The only nonlinearity is the
Fresnel integrals C(u), S(u), which are computed by fixed series / rational
approximations — same per-call cost regardless of inputs.
"""

import math
from scipy.special import fresnel    # for reference; replaceable by polynomial approx
from scipy.integrate import quad


def base_integrals(s, k0, k1, L):
    """Compute (Ic(s), Is(s)) without the theta0 rotation."""
    alpha = (k1 - k0) / L
    if abs(alpha) < 1e-14:
        # Pure arc or straight
        if abs(k0) < 1e-14:
            return s, 0.0
        return math.sin(k0 * s) / k0, (1 - math.cos(k0 * s)) / k0

    sigma = 1.0 if alpha > 0 else -1.0
    a = abs(alpha)
    scale = math.sqrt(math.pi / a)
    gamma = -k0 * k0 / (2 * alpha)

    # u argument bounds
    u_lo = math.sqrt(a / math.pi) * (k0 / alpha)             # at t=0
    u_hi = math.sqrt(a / math.pi) * (s + k0 / alpha)         # at t=s
    # Fresnel: scipy returns (S, C)
    S_hi, C_hi = fresnel(u_hi)
    S_lo, C_lo = fresnel(u_lo)
    dC = C_hi - C_lo
    dS = sigma * (S_hi - S_lo)
    Ic = scale * (math.cos(gamma) * dC - math.sin(gamma) * dS)
    Is = scale * (math.sin(gamma) * dC + math.cos(gamma) * dS)
    return Ic, Is


def clothoid(P0, P1, k0, k1, L):
    """Returns a callable s -> (x,y) and the recovered theta0."""
    A, B = base_integrals(L, k0, k1, L)
    Cx, Cy = P1[0] - P0[0], P1[1] - P0[1]
    # Consistency check
    chord_sq = Cx*Cx + Cy*Cy
    base_sq = A*A + B*B
    if not math.isclose(chord_sq, base_sq, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError(
            f"Inputs inconsistent: |chord|^2={chord_sq:.6g} but A^2+B^2={base_sq:.6g}. "
            "P0,P1,k0,k1,L over-determine the clothoid."
        )
    theta0 = math.atan2(Cy, Cx) - math.atan2(B, A)

    def point(s):
        Ic, Is_ = base_integrals(s, k0, k1, L)
        c, si = math.cos(theta0), math.sin(theta0)
        return (P0[0] + c * Ic - si * Is_,
                P0[1] + si * Ic + c * Is_)
    return point, theta0


def reference_clothoid(P0, theta0, k0, k1, L):
    """Direct numerical-quadrature reference for comparison."""
    alpha = (k1 - k0) / L
    def heading(t): return theta0 + k0*t + 0.5*alpha*t*t
    def point(s):
        cx, _ = quad(lambda t: math.cos(heading(t)), 0, s)
        cy, _ = quad(lambda t: math.sin(heading(t)), 0, s)
        return P0[0] + cx, P0[1] + cy
    return point


# ---------- Verification on several configurations --------------------------
def run_case(name, P0, theta0_true, k0, k1, L):
    ref = reference_clothoid(P0, theta0_true, k0, k1, L)
    P1 = ref(L)
    pts, theta0_rec = clothoid(P0, P1, k0, k1, L)
    # Compare at multiple s values
    max_err = 0.0
    for i in range(11):
        s = L * i / 10
        rx, ry = ref(s)
        ax, ay = pts(s)
        err = math.hypot(rx - ax, ry - ay)
        max_err = max(max_err, err)
    print(f"{name:35s}  theta0 true={theta0_true:+.6f}  rec={theta0_rec:+.6f}  "
          f"max|err|={max_err:.2e}")


print("Deterministic clothoid algorithm — verification\n" + "-"*60)
run_case("clothoid k0=0.0 k1=0.5",         (0, 0), 0.3,   0.00,  0.50, 4.0)
run_case("clothoid k0=-0.3 k1=0.7",        (1, 2), -0.5, -0.30,  0.70, 5.0)
run_case("clothoid k0=0.7 k1=-0.7",        (0, 0), 1.0,   0.70, -0.70, 3.0)
run_case("pure arc k0=k1=0.4 (alpha=0)",   (0, 0), 0.0,   0.40,  0.40, 2.5)
run_case("straight k0=k1=0",               (0, 0), 1.2,   0.00,  0.00, 7.0)
run_case("near-straight k0=1e-9 k1=2e-9",  (0, 0), 0.1,   1e-9,  2e-9, 1.0)
run_case("large curvature k0=2 k1=-1",     (3,-4), -2.0,  2.0,  -1.0, 1.5)
