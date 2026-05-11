# SPDX-FileCopyrightText: 2026 Merkator Group
# SPDX-License-Identifier: LicenseRef-Merkator-Proprietary-NoAITraining
"""
Run the Halley / Newton Python solvers over the full golden-vector
corpus and emit the same JSON shape as the C# / Java / TypeScript
benchmarks. Intended to be consumed by run_all_benches.py.
"""
from __future__ import annotations

import json
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from clothoid_halleyL_bench import solve_halley_L, solve_newton_L  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
GOLDEN = ROOT / "data" / "golden_vectors.json"

WARMUPS = 5
REPEATS = 50


def _run(cases, solver):
    iters = 0
    for c in cases:
        _, it = solver(c["P0"], c["P1"], c["k0"], c["k1"])
        iters += it
    return iters


def main() -> int:
    with GOLDEN.open(encoding="utf-8") as f:
        cases = json.load(f)["cases"]

    for _ in range(WARMUPS):
        _run(cases, solve_halley_L)
        _run(cases, solve_newton_L)

    halley_ms = []
    newton_ms = []
    halley_iters = newton_iters = 0
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        halley_iters = _run(cases, solve_halley_L)
        halley_ms.append((time.perf_counter() - t0) * 1000.0)

        t0 = time.perf_counter()
        newton_iters = _run(cases, solve_newton_L)
        newton_ms.append((time.perf_counter() - t0) * 1000.0)

    halley_ms.sort()
    newton_ms.sort()
    median_h = halley_ms[REPEATS // 2]
    median_n = newton_ms[REPEATS // 2]
    min_h    = halley_ms[0]
    min_n    = newton_ms[0]

    result = {
        "language":         "Python 3 (NumPy)",
        "runtime":          f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "cases":            len(cases),
        "halley_us":        1000.0 * median_h / len(cases),
        "newton_us":        1000.0 * median_n / len(cases),
        "halley_us_min":    1000.0 * min_h    / len(cases),
        "newton_us_min":    1000.0 * min_n    / len(cases),
        "halley_iter_mean": halley_iters / len(cases),
        "newton_iter_mean": newton_iters / len(cases),
    }
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
