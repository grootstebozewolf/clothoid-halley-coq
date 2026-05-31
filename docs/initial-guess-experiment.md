# Initial-guess sensitivity and the dimensionless convergence field

Design note for the experiment hardening Section 5 of
*A Verified-Azimuth, Zero-Friction Halley Solver for the Chord-Length
Parameter L*. Harness: [`python/clothoid_initial_guess_experiment.py`](../python/clothoid_initial_guess_experiment.py).

## Motivation

The referee's central observation is that iteration count is dominated by
initial-guess quality far more than by the Halley-vs-Newton order, and that the
reported mean of 2.28 may "describe the dataset as much as the method." This
experiment settles that quantitatively by exploiting the symmetries of the
residual, so the entire question reduces to a single 2-D field that can be
mapped exhaustively.

## Symmetry stack

The residual is `f(L) = L²(P²+Q²) − d²` with the six moments built from
`ψ(τ) = κ₀τ + (κ₁−κ₀)τ²/2`. Three nested symmetries collapse the parameter
space:

1. **SE(2) invariance.** `f` depends on the geometry only through the chord
   length `d` and the curvatures `(κ₀,κ₁)`; position and orientation drop out.
   (5 → 3 parameters.)

2. **Scaling.** `L·ψ(τ) = aτ + (b−a)τ²/2` with the **dimensionless curvatures**
   `a = κ₀L`, `b = κ₁L`. Hence `r = √(P²+Q²)` is a function `r(a,b)` and
   `d/L = r(a,b)`. The iteration count is therefore a scalar field `N(a*,b*)`
   over the monotone square `[−π,π]²`, where `(a*,b*) = (κ₀L*, κ₁L*)` is the
   dimensionless solution. (3 → 2.)

3. **Klein-four group.** `r(a,b) = r(b,a) = r(−a,−b) = r(−b,−a)`
   (arc reversal `τ→1−τ`, and the mirror clothoid `Q→−Q`). `N` inherits the
   group, so only the wedge `{a ≥ |b|}` — one quarter of the square — need be
   computed. (4× reduction.)

Net: 5 physical parameters → a 2-D field on ¼ of `[−π,π]²`.

## Hypotheses (all confirmed)

| # | Hypothesis | Result |
|---|------------|--------|
| H1 | `N` depends only on `(a*,b*)`, not on `(κ₀,κ₁,L*)` separately | Holds for `d²≫1`; broken by at most ±1 iteration below ~1 m by the `max(d²,1)` tolerance floor |
| H2 | `N(a,b)=N(b,a)=N(−a,−b)=N(−b,−a)` | Confirmed; symmetry-fill vs direct recompute differs by ≤1 iteration on a thin threshold set (machine-precision quadrature asymmetry tipping the stopping test), **3.97× compute saving** |
| H3 | `N` collapses onto the initial relative error `1−r(a,b)` | Confirmed (see `initguess_collapse.png`) |
| H4 | soundness/large-`N` confined to a boundary layer near `|a|,|b|=π` | Confirmed; interior is benign, `N_max = 14` (Halley), `21` (Newton) at the corners |

## Key findings

- **The ProRail corpus is near-straight in dimensionless terms.** The full
  9,058-record corpus occupies `|a| ≤ 0.38`, `|b| ≤ 0.43` — about 12% of the
  way to the `π` boundary. The harness independently reproduces the published
  corpus means (Halley **2.28**, Newton **2.92**; per-case agreement within 1
  iteration, the ±1 from RD-coordinate cancellation in `d`).
- **`L₀ = d` is near-optimal, and the low count is the initial guess.** Over a
  uniform sweep of the square the mean iteration count is *minimised* at the
  policy `L₀ = c·d` with `c = 1` (3.46 iterations). On the corpus the mean rises
  from **2.28** at `c=1` to **4.0** at `c=0.75` or `c=1.25` — a ±25% nudge of the
  initial guess nearly doubles the iteration count. The 2.28 figure is thus a
  property of `L₀=d` landing almost on the root in the near-straight railway
  regime, not a generic property of Halley.

## Reproducing

```
pip install numpy scipy matplotlib
python3 python/clothoid_initial_guess_experiment.py
```

Drives the real solver (`moments`, `solve_halley_L`, `solve_newton_L` from
`clothoid_halleyL_bench`); an `L₀`-parameterised variant is asserted to
reproduce the shipped solver exactly when `L₀=d`. Outputs land in
`docs/mathematics/`: `initguess_iterations.png`, `initguess_collapse.png`,
`initguess_l0_sensitivity.png`, `initguess_scaling.png`, and the machine-readable
`initguess_results.json`.
