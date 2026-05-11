"""
Drive the four language benchmarks (Python / C# / Java / TypeScript)
sequentially, collect their JSON output, and write a unified
benchmark_results.json that the paper / bar-chart generator consume.

Each language bench prints exactly one JSON object on stdout matching:

  {
    "language":         "...",
    "runtime":          "...",
    "cases":            int,
    "halley_us":        float,    # median per-case microseconds
    "newton_us":        float,
    "halley_us_min":    float,
    "newton_us_min":    float,
    "halley_iter_mean": float,
    "newton_iter_mean": float
  }
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "data" / "benchmark_results.json"


# (label, cwd, command list) tuples.
BENCHES = [
    ("python", ROOT / "python",
     [sys.executable, "bench_corpus.py"]),
    ("csharp", ROOT / "csharp" / "Clothoid.Halley.Bench",
     ["dotnet", "run", "-c", "Release", "-nologo", "--verbosity", "quiet"]),
    ("java",   ROOT,
     ["java", "-jar", str(ROOT / "java" / "target" / "clothoid-halley-1.0.0.jar")]),
    ("typescript", ROOT / "typescript",
     ["npm", "run", "bench", "--silent"]),
]


def _parse_json_block(text: str) -> dict:
    """Pull the first {...} JSON object out of a stdout transcript."""
    m = re.search(r"\{[^{}]*?(?:\{[^{}]*\}[^{}]*?)*\}", text, re.DOTALL)
    if not m:
        raise ValueError("no JSON block found in stdout")
    return json.loads(m.group(0))


def main() -> int:
    results = []
    for label, cwd, cmd in BENCHES:
        print(f"[{label}] running {' '.join(cmd)} ...", file=sys.stderr)
        t0 = time.time()
        # Windows: npm/dotnet are .cmd shims, need shell=True to resolve them.
        use_shell = sys.platform == "win32"
        cmd_arg = " ".join(f'"{c}"' if " " in c else c for c in cmd) if use_shell else cmd
        proc = subprocess.run(cmd_arg, cwd=cwd, capture_output=True, text=True, shell=use_shell)
        dt = time.time() - t0
        if proc.returncode != 0:
            sys.stderr.write(f"[{label}] FAILED exit={proc.returncode}\n--- stderr ---\n{proc.stderr}\n")
            return 1
        try:
            obj = _parse_json_block(proc.stdout)
        except Exception as e:
            sys.stderr.write(f"[{label}] could not parse JSON: {e}\n--- stdout ---\n{proc.stdout}\n")
            return 1
        obj["wall_s"] = dt
        results.append(obj)
        print(f"[{label}] Halley {obj['halley_us']:.2f} us/op  Newton {obj['newton_us']:.2f} us/op  "
              f"({dt:.1f}s wall)", file=sys.stderr)

    payload = {
        "generated":      time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "platform":       sys.platform,
        "corpus_source":  "data/golden_vectors.json (real ProRail Spoorgeometrie clothoids)",
        "corpus_size":    results[0]["cases"],
        "results":        results,
    }
    with OUT.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(f"\nwrote {OUT.relative_to(ROOT)}", file=sys.stderr)

    # Print a console-friendly table
    print("\n  language               Halley us/op   Newton us/op   Halley iters   Newton iters")
    print(  "  " + "-" * 84)
    for r in results:
        print(f"  {r['language']:22s} {r['halley_us']:11.2f}   {r['newton_us']:11.2f}   "
              f"{r['halley_iter_mean']:11.2f}   {r['newton_iter_mean']:11.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
