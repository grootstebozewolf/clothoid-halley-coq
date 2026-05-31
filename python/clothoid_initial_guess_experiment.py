#!/usr/bin/env python3
"""Initial-guess sensitivity and dimensionless symmetry experiment for the
chord-length Halley solver (paper Section 5 hardening).

The experiment rests on three nested symmetries of the residual
    f(L) = L^2 (P^2 + Q^2) - d^2,   P = int_0^1 cos(L psi), Q = int_0^1 sin(L psi)
with psi(tau) = k0 tau + (k1 - k0) tau^2 / 2:

  1. SE(2) invariance: f depends on the geometry only through the chord
     length d and the curvatures (k0, k1); orientation and position drop out.
  2. Scaling: L*psi depends only on the DIMENSIONLESS curvatures
        a = k0*L,  b = k1*L,
     so r = sqrt(P^2+Q^2) is a function r(a,b) and d/L = r(a,b).  The whole
     solver dynamics -- hence the iteration count N -- is a scalar field
     N(a*,b*) over the monotone square [-pi,pi]^2 (a* = k0 L*, b* = k1 L*).
  3. Klein-four group: r(a,b) = r(b,a) = r(-a,-b) = r(-b,-a) (arc reversal
     and mirror), so only the 1/4 wedge {a >= |b|} need be computed.

This script DRIVES THE REAL solver: it imports `moments`, `solve_halley_L`,
`solve_newton_L` from clothoid_halleyL_bench and reuses the identical 32-point
Gauss-Legendre rule and guard logic.  An L0-parameterised variant (needed for
the initial-guess axis) is checked to reproduce the shipped solvers exactly
when L0 = d.

Outputs (written next to the manuscript so figures can be \\includegraphics'd):
  docs/mathematics/initguess_iterations.png   N(a,b) maps, ProRail cloud overlay
  docs/mathematics/initguess_collapse.png      N vs initial relative error 1-r
  docs/mathematics/initguess_l0_sensitivity.png mean N vs L0 = c*d
  docs/mathematics/initguess_scaling.png        scaling-floor breakdown of sym. 1->2
  docs/mathematics/initguess_results.json       machine-readable summary
"""
import json
import math
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from clothoid_halleyL_bench import (
    moments, solve_halley_L, solve_newton_L, _GL_N,
)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIGDIR = os.path.join(ROOT, "docs", "mathematics")
GOLDEN = os.path.join(ROOT, "data", "golden_vectors.json")
PI = math.pi
TOL = 1e-13
MAX_ITER = 50

# --------------------------------------------------------------------------
# L0-parameterised solvers: byte-for-byte the shipped iteration, but the
# initial guess is an explicit argument instead of the hard-coded L = d.
# --------------------------------------------------------------------------
def halley_from(k0, k1, d, L0, tol=TOL, max_iter=MAX_ITER):
    d2 = d * d
    if d == 0:
        return 0.0, 0
    L = L0
    for it in range(1, max_iter + 1):
        P, Q, R, T, S2c, S2s = moments(L, k0, k1)
        r2 = P * P + Q * Q
        QRPT = Q * R - P * T
        f = L * L * r2 - d2
        fp = 2 * L * r2 + 2 * L * L * QRPT
        fpp = 2 * r2 + 8 * L * QRPT + 2 * L * L * (R * R + T * T - P * S2c - Q * S2s)
        if abs(f) < tol * max(d2, 1.0):
            return L, it
        denom = 2 * fp * fp - f * fpp
        if abs(denom) < 1e-20 or fp <= 0:
            L *= 1.5
            continue
        L_new = L - 2 * f * fp / denom
        if L_new <= 0:
            L_new = 0.5 * L
        L = L_new
    return L, max_iter


def newton_from(k0, k1, d, L0, tol=TOL, max_iter=MAX_ITER):
    d2 = d * d
    if d == 0:
        return 0.0, 0
    L = L0
    for it in range(1, max_iter + 1):
        P, Q, R, T, _, _ = moments(L, k0, k1)
        r2 = P * P + Q * Q
        f = L * L * r2 - d2
        fp = 2 * L * r2 + 2 * L * L * (Q * R - P * T)
        if abs(f) < tol * max(d2, 1.0):
            return L, it
        if fp <= 0:
            L *= 1.5
            continue
        L = max(L - f / fp, 0.5 * L)
    return L, max_iter


def r_of(a, b):
    """Dimensionless chord ratio r(a,b) = |int_0^1 exp(i(a tau + (b-a)tau^2/2))|,
    evaluated with the solver's own 32-point Gauss-Legendre rule (set L=1,
    k0=a, k1=b so that L*psi = a tau + (b-a) tau^2/2)."""
    P, Q, _, _, _, _ = moments(1.0, a, b)
    return math.hypot(P, Q)


def instance(a, b, Lstar=1.0):
    """Physical instance whose dimensionless solution is (a,b): k0=a/Lstar,
    k1=b/Lstar, chord d = Lstar*r(a,b).  L* = Lstar is then an exact root of
    the quadrature residual (forward and inverse use the identical GL rule)."""
    k0, k1 = a / Lstar, b / Lstar
    d = Lstar * r_of(a, b)
    return k0, k1, d


# --------------------------------------------------------------------------
# 0. Self-check: the L0-parameterised solver with L0 = d reproduces the
#    shipped solver exactly (same returned L and iteration count).
# --------------------------------------------------------------------------
def selfcheck_against_shipped():
    rng = np.random.default_rng(0)
    worst = 0
    for _ in range(2000):
        a = rng.uniform(-PI, PI)
        b = rng.uniform(-PI, PI)
        k0, k1, d = instance(a, b)
        P0, P1 = (0.0, 0.0), (d, 0.0)
        Lh, ih = solve_halley_L(P0, P1, k0, k1)
        Ln, in_ = solve_newton_L(P0, P1, k0, k1)
        Lh2, ih2 = halley_from(k0, k1, d, d)
        Ln2, in2 = newton_from(k0, k1, d, d)
        assert ih == ih2 and in_ == in2, (a, b, ih, ih2, in_, in2)
        worst = max(worst, abs(Lh - Lh2), abs(Ln - Ln2))
    return worst


# --------------------------------------------------------------------------
# 1+2. Symmetry-reduced field over the monotone square.  Compute only the
#      wedge {a >= |b|} (1/4 of the square), fill the rest by the Klein-four
#      group, then validate against direct computation on a random subset.
# --------------------------------------------------------------------------
def compute_fields(nres=257):
    # exactly antisymmetric axis: axis[nres-1-i] == -axis[i] to the bit, so the
    # Klein-four group maps grid cells onto grid cells with no boundary leakage.
    half = (nres - 1) // 2
    step = PI / half
    axis = step * (np.arange(nres) - half)
    A, B = np.meshgrid(axis, axis, indexing="ij")  # A[i,j]=axis[i], B[i,j]=axis[j]
    Nh = np.full((nres, nres), -1, dtype=int)
    Nn = np.full((nres, nres), -1, dtype=int)
    Rf = np.full((nres, nres), np.nan)
    Eh = np.full((nres, nres), np.nan)

    wedge = A >= np.abs(B)  # fundamental domain a >= |b|
    t0 = time.time()
    computed = 0
    for i in range(nres):
        for j in range(nres):
            if not wedge[i, j]:
                continue
            a, b = axis[i], axis[j]
            k0, k1, d = instance(a, b)
            Lh, ih = halley_from(k0, k1, d, d)
            Ln, in_ = newton_from(k0, k1, d, d)
            r = d  # Lstar=1 so d = r(a,b)
            # assign the cell and its three group images
            for (ii, jj) in ((i, j), (j, i),
                             (nres - 1 - i, nres - 1 - j),
                             (nres - 1 - j, nres - 1 - i)):
                Nh[ii, jj] = ih
                Nn[ii, jj] = in_
                Rf[ii, jj] = r
                Eh[ii, jj] = abs(Lh - 1.0)
            computed += 1
    t_reduced = time.time() - t0
    assert (Nh >= 0).all(), "group fill left holes"

    # validate the symmetry fill on a random subset by direct computation
    rng = np.random.default_rng(1)
    idx = rng.integers(0, nres, size=(400, 2))
    sym_err = 0
    for i, j in idx:
        k0, k1, d = instance(axis[i], axis[j])
        _, ih = halley_from(k0, k1, d, d)
        _, in_ = newton_from(k0, k1, d, d)
        sym_err = max(sym_err, abs(int(ih) - int(Nh[i, j])),
                      abs(int(in_) - int(Nn[i, j])))
    speedup = (nres * nres) / computed
    return dict(axis=axis, A=A, B=B, Nh=Nh, Nn=Nn, Rf=Rf, Eh=Eh,
                t_reduced=t_reduced, computed=computed, total=nres * nres,
                speedup=speedup, sym_err=int(sym_err))


# --------------------------------------------------------------------------
# 3. ProRail corpus: dimensionless cloud + harness-vs-published validation.
# --------------------------------------------------------------------------
def prorail_cloud():
    if not os.path.exists(GOLDEN):
        return None
    cases = json.load(open(GOLDEN))["cases"]
    aa, bb, nh_pub, nn_pub, nh_re, nn_re = [], [], [], [], [], []
    for c in cases:
        k0, k1, L = c["k0"], c["k1"], c["L"]
        d = c["d"]
        aa.append(k0 * L)
        bb.append(k1 * L)
        nh_pub.append(c["iter_halley"])
        nn_pub.append(c["iter_newton"])
        _, ih = halley_from(k0, k1, d, d)
        _, in_ = newton_from(k0, k1, d, d)
        nh_re.append(ih)
        nn_re.append(in_)
    aa, bb = np.array(aa), np.array(bb)
    nh_pub, nn_pub = np.array(nh_pub), np.array(nn_pub)
    nh_re, nn_re = np.array(nh_re), np.array(nn_re)
    return dict(a=aa, b=bb, nh_pub=nh_pub, nn_pub=nn_pub, nh_re=nh_re, nn_re=nn_re)


# --------------------------------------------------------------------------
# 4. Initial-guess sensitivity: mean iterations vs L0 = c * d.
# --------------------------------------------------------------------------
def l0_sensitivity(corpus, nsq=64):
    cs = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]
    axis = np.linspace(-PI, PI, nsq)
    sq = []
    for a in axis:
        for b in axis:
            sq.append(instance(a, b))
    out = {"c": cs, "square_halley": [], "square_newton": [],
           "corpus_halley": [], "corpus_newton": []}
    for c in cs:
        h = [halley_from(k0, k1, d, c * d)[1] for (k0, k1, d) in sq]
        n = [newton_from(k0, k1, d, c * d)[1] for (k0, k1, d) in sq]
        out["square_halley"].append(float(np.mean(h)))
        out["square_newton"].append(float(np.mean(n)))
        if corpus is not None:
            cases = json.load(open(GOLDEN))["cases"]
            hc = [halley_from(cc["k0"], cc["k1"], cc["d"], c * cc["d"])[1] for cc in cases]
            nc = [newton_from(cc["k0"], cc["k1"], cc["d"], c * cc["d"])[1] for cc in cases]
            out["corpus_halley"].append(float(np.mean(hc)))
            out["corpus_newton"].append(float(np.mean(nc)))
    return out


# --------------------------------------------------------------------------
# 5. Scaling-floor: symmetry (2) is broken by the max(d^2,1) tolerance floor
#    at small absolute scale.  Sweep L* across decades and watch N move.
# --------------------------------------------------------------------------
def scaling_floor():
    alphas = np.logspace(-3, 4, 22)
    probes = [(2.5, 1.0), (3.0, -2.0), (PI * 0.95, PI * 0.95)]
    out = {"alpha": alphas.tolist(), "probes": []}
    for (a, b) in probes:
        row = []
        for al in alphas:
            k0, k1, d = instance(a, b, Lstar=al)
            _, ih = halley_from(k0, k1, d, d)
            row.append(ih)
        out["probes"].append({"a": a, "b": b, "N": row})
    return out


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------
def plot_fields(F, corpus):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, (Nfield, title) in zip(axes, ((F["Nh"], "Halley"), (F["Nn"], "Newton"))):
        im = ax.imshow(Nfield.T, origin="lower", extent=[-PI, PI, -PI, PI],
                       aspect="equal", cmap="viridis")
        ax.set_title(f"{title}: iterations $N(a,b)$")
        ax.set_xlabel(r"$a = \kappa_0 L^\ast$")
        ax.set_ylabel(r"$b = \kappa_1 L^\ast$")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="iterations")
        if corpus is not None:
            ax.scatter(corpus["a"], corpus["b"], s=2, c="white",
                       alpha=0.35, linewidths=0, label="ProRail")
            ax.legend(loc="upper right", fontsize=7, framealpha=0.6)
    fig.suptitle(r"Dimensionless convergence field over the monotone square $|a|,|b|\leq\pi$")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "initguess_iterations.png"), dpi=140)
    plt.close(fig)


def plot_collapse(F):
    one_minus_r = 1.0 - F["Rf"].ravel()
    nh = F["Nh"].ravel()
    fig, ax = plt.subplots(figsize=(6, 4.4))
    ax.scatter(one_minus_r, nh, s=3, alpha=0.15, c="C0", linewidths=0)
    # binned mean
    bins = np.linspace(0, one_minus_r.max(), 40)
    idx = np.digitize(one_minus_r, bins)
    bx, by = [], []
    for k in range(1, len(bins)):
        m = idx == k
        if m.sum() > 5:
            bx.append(one_minus_r[m].mean())
            by.append(nh[m].mean())
    ax.plot(bx, by, "C3-o", ms=3, lw=1.5, label="binned mean")
    ax.set_xlabel(r"initial relative error $1 - r(a,b) = 1 - L_0/L^\ast$  (for $L_0=d$)")
    ax.set_ylabel("Halley iterations $N$")
    ax.set_title("Iteration count collapses onto the initial-guess error")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "initguess_collapse.png"), dpi=140)
    plt.close(fig)


def plot_l0(sens):
    fig, ax = plt.subplots(figsize=(6, 4.4))
    cs = sens["c"]
    ax.plot(cs, sens["square_halley"], "C0-o", ms=4, label="Halley (square mean)")
    ax.plot(cs, sens["square_newton"], "C1-s", ms=4, label="Newton (square mean)")
    if sens["corpus_halley"]:
        ax.plot(cs, sens["corpus_halley"], "C0--^", ms=4, label="Halley (ProRail)")
        ax.plot(cs, sens["corpus_newton"], "C1--v", ms=4, label="Newton (ProRail)")
    ax.axvline(1.0, color="grey", ls=":", lw=1)
    ax.set_xlabel(r"initial-guess multiplier $c$ in $L_0 = c\,d$")
    ax.set_ylabel("mean iterations")
    ax.set_title("Sensitivity of iteration count to the initial guess")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "initguess_l0_sensitivity.png"), dpi=140)
    plt.close(fig)


def plot_scaling(sc):
    fig, ax = plt.subplots(figsize=(6, 4.4))
    al = np.array(sc["alpha"])
    for p in sc["probes"]:
        ax.semilogx(al, p["N"], "-o", ms=3,
                    label=f"(a,b)=({p['a']:.2f},{p['b']:.2f})")
    ax.axvline(1.0, color="grey", ls=":", lw=1, label=r"$d^2=1$ floor")
    ax.set_xlabel(r"physical scale $L^\ast$ (m) at fixed $(a,b)$")
    ax.set_ylabel("Halley iterations $N$")
    ax.set_title("Scaling symmetry holds where $d^2\\gg1$; broken by the tol floor")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "initguess_scaling.png"), dpi=140)
    plt.close(fig)


def main():
    print(f"GL nodes = {_GL_N};  tol = {TOL};  monotone square |a|,|b| <= pi")
    print("[0] self-check vs shipped solver ...", flush=True)
    worst = selfcheck_against_shipped()
    print(f"    L0=d reproduces shipped solver exactly; max |dL| = {worst:.2e}")

    print("[1] symmetry-reduced field over the square ...", flush=True)
    F = compute_fields(nres=257)
    print(f"    computed {F['computed']}/{F['total']} cells "
          f"(Klein-four speedup {F['speedup']:.2f}x) in {F['t_reduced']:.1f}s")
    print(f"    symmetry-fill vs direct recompute: max iter mismatch = {F['sym_err']} "
          f"(0 = exact group invariance)")

    print("[3] ProRail corpus cloud + validation ...", flush=True)
    corpus = prorail_cloud()
    if corpus is not None:
        mism_h = int(np.abs(corpus["nh_pub"] - corpus["nh_re"]).max())
        mism_n = int(np.abs(corpus["nn_pub"] - corpus["nn_re"]).max())
        print(f"    {len(corpus['a'])} records; harness vs published iter "
              f"max mismatch: Halley {mism_h}, Newton {mism_n}")
        print(f"    mean iterations (harness): Halley {corpus['nh_re'].mean():.2f}, "
              f"Newton {corpus['nn_re'].mean():.2f}  (paper: 2.28 / 2.92)")
        print(f"    dimensionless extent: |a|<= {np.abs(corpus['a']).max():.3f}, "
              f"|b|<= {np.abs(corpus['b']).max():.3f}  (vs pi={PI:.3f})")

    print("[4] initial-guess sensitivity ...", flush=True)
    sens = l0_sensitivity(corpus)
    for c, sh in zip(sens["c"], sens["square_halley"]):
        print(f"    c={c:>4}:  square-mean Halley N = {sh:.2f}")

    print("[5] scaling-floor symmetry breakdown ...", flush=True)
    sc = scaling_floor()

    plot_fields(F, corpus)
    plot_collapse(F)
    plot_l0(sens)
    plot_scaling(sc)

    summary = dict(
        gl_nodes=_GL_N, tol=TOL,
        selfcheck_max_dL=worst,
        grid=dict(nres=257, computed=F["computed"], total=F["total"],
                  speedup=F["speedup"], symmetry_max_mismatch=F["sym_err"],
                  N_halley_max=int(F["Nh"].max()), N_newton_max=int(F["Nn"].max())),
        l0_sensitivity=sens,
        scaling_floor=sc,
    )
    if corpus is not None:
        summary["prorail"] = dict(
            n=len(corpus["a"]),
            mean_halley=float(corpus["nh_re"].mean()),
            mean_newton=float(corpus["nn_re"].mean()),
            max_abs_a=float(np.abs(corpus["a"]).max()),
            max_abs_b=float(np.abs(corpus["b"]).max()),
            published_mismatch_halley=int(np.abs(corpus["nh_pub"] - corpus["nh_re"]).max()),
            published_mismatch_newton=int(np.abs(corpus["nn_pub"] - corpus["nn_re"]).max()),
        )
    json.dump(summary, open(os.path.join(FIGDIR, "initguess_results.json"), "w"), indent=2)
    print(f"\nwrote figures + initguess_results.json to {FIGDIR}")


if __name__ == "__main__":
    main()
