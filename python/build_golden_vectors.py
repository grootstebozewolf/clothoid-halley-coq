"""
Build the golden-vector test suite that every implementation
(Python / C# / Java / TypeScript) must reproduce numerically.

Reads  data/prorail_clothoids.json.gz       (raw ProRail snapshot)
Writes data/golden_vectors.json             (filtered + solved corpus)

Filtering rules (matching the monotone-branch assumption in the paper):
- chord d > 0
- |k0 * L_design| <= pi
- |k1 * L_design| <= pi

For each surviving record, both solve_halley_L and solve_newton_L are
run; the record is recorded if both converged inside max_iter and
agreed within 1e-9 metres.

Output JSON shape:

    {
      "source":          "data/prorail_clothoids.json.gz",
      "source_fetched":  "<UTC timestamp>",
      "tol":             1e-13,
      "max_iter":        50,
      "count":           <int>,
      "cases": [
        { "objectid": int,
          "P0":       [x, y],
          "P1":       [x, y],
          "k0":       float,
          "k1":       float,
          "L_design": float,
          "d":        float,
          "L":        float,   # solver-agreed L (Halley, == Newton within 1e-9)
          "iter_halley": int,
          "iter_newton": int
        }, ...
      ]
    }
"""
from __future__ import annotations

import gzip
import json
import math
import sys
from pathlib import Path

# import the solver from the same package
sys.path.insert(0, str(Path(__file__).resolve().parent))
from clothoid_halleyL_bench import solve_halley_L, solve_newton_L  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
IN   = ROOT / "data" / "prorail_clothoids.json.gz"
OUT  = ROOT / "data" / "golden_vectors.json"

TOL_CONVERGE = 1e-13
MAX_ITER     = 50
AGREE_TOL    = 1e-9


def main() -> int:
    with gzip.open(IN, "rt", encoding="utf-8") as g:
        raw = json.load(g)
    records = raw["records"]
    sys.stderr.write(f"loaded {len(records)} raw records from {IN.name}\n")

    cases:    list[dict] = []
    skipped = {
        "zero_chord":        0,
        "monotone_violation": 0,
        "halley_no_converge": 0,
        "newton_no_converge": 0,
        "method_disagree":    0,
    }

    for r in records:
        P0 = r["P0"]; P1 = r["P1"]
        k0 = r["k0"]; k1 = r["k1"]
        L_des = r["L_design"]
        Cx, Cy = P1[0] - P0[0], P1[1] - P0[1]
        d = math.hypot(Cx, Cy)
        if d <= 0:
            skipped["zero_chord"] += 1
            continue
        if abs(k0 * L_des) > math.pi or abs(k1 * L_des) > math.pi:
            skipped["monotone_violation"] += 1
            continue
        LH, iH = solve_halley_L(P0, P1, k0, k1, tol=TOL_CONVERGE, max_iter=MAX_ITER)
        LN, iN = solve_newton_L(P0, P1, k0, k1, tol=TOL_CONVERGE, max_iter=MAX_ITER)
        if iH >= MAX_ITER:
            skipped["halley_no_converge"] += 1
            continue
        if iN >= MAX_ITER:
            skipped["newton_no_converge"] += 1
            continue
        if abs(LH - LN) > AGREE_TOL:
            skipped["method_disagree"] += 1
            continue
        cases.append({
            "objectid":     r["objectid"],
            "P0":           P0,
            "P1":           P1,
            "k0":           k0,
            "k1":           k1,
            "L_design":     L_des,
            "d":            d,
            "L":            LH,
            "iter_halley":  iH,
            "iter_newton":  iN,
        })

    sys.stderr.write("\nfilter / solve summary:\n")
    for k, v in skipped.items():
        sys.stderr.write(f"  {k:20s} {v:6d}\n")
    sys.stderr.write(f"  {'usable cases':20s} {len(cases):6d}\n")

    payload = {
        "source":         f"data/{IN.name}",
        "source_fetched": raw.get("fetched"),
        "tol":            TOL_CONVERGE,
        "max_iter":       MAX_ITER,
        "count":          len(cases),
        "cases":          cases,
    }
    with OUT.open("w", encoding="utf-8") as f:
        json.dump(payload, f, separators=(",", ":"))
    sys.stderr.write(f"\nwrote {OUT.relative_to(ROOT)}  ({OUT.stat().st_size/1024:.1f} KB)\n")

    # quick aggregate stats
    if cases:
        import statistics
        iH = [c["iter_halley"] for c in cases]
        iN = [c["iter_newton"] for c in cases]
        sys.stderr.write("\niteration counts (Halley / Newton):\n")
        sys.stderr.write(f"  mean    {statistics.mean(iH):6.3f}  /  {statistics.mean(iN):6.3f}\n")
        sys.stderr.write(f"  median  {statistics.median(iH):6.0f}  /  {statistics.median(iN):6.0f}\n")
        sys.stderr.write(f"  max     {max(iH):6d}  /  {max(iN):6d}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
