# clothoid-halley-coq

Halley iteration for clothoid (Euler-spiral) fitting, with a Coq formalization
of the key derivative identities under Coquelicot.

## What's here

### `coq/`

Three Coq files, all compile cleanly under Coq 8.13.1 + Coquelicot 3.x and pass
`coqchk` with no `type-in-type`, no unsafe (co)fixpoints, no positivity holes.

- **`Clothoid.v`** — Bertolazzi–Frego G¹ residual `f(A) = Y₀(2A, δ−A, φ₀)`.
  Proves `f'(A) = ∫(t²−t)·cos(φ)` and `f''(A) = −∫(t²−t)²·sin(φ)` via
  parameter-differentiation under the integral (Coquelicot's
  `is_derive_RInt_param_aux`).
- **`Clothoid_L.v`** — chord-length residual `f(L) = L²·(P²+Q²) − d²` viewed
  as a function of L. Proves the four integral-level derivatives
  `P'=−T`, `Q'=R`, `R'=−S2s`, `T'=S2c` and the first-derivative composite
  `f'(L) = 2L(P²+Q²) + 2L²(QR−PT)`. One lemma (`f_L_second_eq`) remains
  admitted; details below.
- **`ClothoidPolish.v`** — moment-notation polish: defines `f`, `Xₖ`, `Yₖ`
  per Bertolazzi–Frego and proves `f'(A) = X₂−X₁` and the integral identity
  `−(Y₄−2Y₃+Y₂) = ∫−(t²−t)²·sin(φ)` (`Y_pattern_eq`).

#### Build

```
coqc Clothoid.v
coqc Clothoid_L.v
coqc -R . "" ClothoidPolish.v
coqchk -silent -o Clothoid Clothoid_L
coqchk -silent -o -R . "" ClothoidPolish
```

#### Axioms

All proved theorems use only Coq's standard four axioms (classical logic,
decidable Dedekind reals, functional extensionality) — the same axioms
underlying Coquelicot's analysis library. Verified with `Print Assumptions`.

#### The one remaining admit

`f_L_second_eq` in `Clothoid_L.v` — the composite second derivative on L. The
analytical content (the four integral-level derivatives `is_derive_P/Q/R/T`)
is fully proved; the admit is mechanical product-rule composition. The
obstacle is Coq's typeclass-vs-Reals unifier hitting deeper-than-first-order
nestings: Coquelicot's `is_derive_mult` etc. use generic `plus`/`mult` while
expressions use `Rplus`/`Rmult`. The first-derivative composite escapes this
by careful explicit instances; the second derivative's deeper tree does not.
A clean fix is documented in the source: ~100–200 lines of explicit
instances or a dedicated `Ltac` for derive-composition.

### `python/`

Numerical experiments and benchmarks behind the algorithmic claims.

- **Sympy derivations**
  - `clothoid_derivation.py` — closed-form `θ₀ = atan2(C_y,C_x) − atan2(B,A)`
    when `(P₀,P₁,κ₀,κ₁,L)` are given.
  - `clothoid_newton_derivation.py` — Newton on the chord-length residual.
  - `clothoid_halley_derive.py` — sympy-verified f', f'' for the
    Bertolazzi–Frego residual.
  - `clothoid_halleyL_derive.py` — sympy-verified f', f'' for the
    chord-length residual on L.

- **Verifications and benchmarks**
  - `clothoid_verify.py` — closed-form algorithm vs reference quadrature.
  - `clothoid_newton_verify.py` — Newton-on-L iteration counts.
  - `clothoid_halley_benchmark.py` — Halley vs Newton on the Bertolazzi–Frego
    kernel, 21×21 (φ₀,φ₁) hypercube.
  - `clothoid_stress.py` — same kernel, 101×101 grid (10 200 samples).
    Result: Halley caps worst case at 3 iterations vs Newton's 4 (~48% of
    samples reduce by one iteration).
  - `clothoid_halleyL_bench.py` — Halley vs Newton on the chord-length-on-L
    kernel, with the *verified* `f''(L)` formula.
  - `peer_review_halley.py` — independent review of a third-party Halley
    implementation; identifies the labelling/derivative defects.

#### Notes on the python work

The Python derivations include the corrected `f''(L)` formula that an earlier
sketch had wrong by three sign/coupling errors. Both the wrong formula and
the corrected one are shown side-by-side in `clothoid_halleyL_derive.py` with
a finite-difference cross-check that rejects the wrong one and accepts the
correct one.

## Status

Substantive work, single-author, private. Not intended for external
distribution in its current form — the publication line in the README of the
worktree describes what would be ready to claim.

## Build environment

- Coq 8.13.1 (February 2021), Coquelicot 3.x, installed at `C:\Coq`.
- Python 3.13/3.14 with numpy, scipy, sympy.
