"""
Render the two benchmark bar-charts that appear in the paper from the
measured numbers in data/benchmark_results.json (produced by
python/run_all_benches.py).

Reads:   data/benchmark_results.json
Writes:  docs/mathematics/iterations_comparison.png
         docs/mathematics/time_comparison.png
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT    = Path(__file__).resolve().parent.parent.parent
INPUT   = ROOT / "data" / "benchmark_results.json"
OUT_DIR = ROOT / "docs" / "mathematics"


def main() -> None:
    with INPUT.open(encoding="utf-8") as f:
        bench = json.load(f)
    results = bench["results"]
    languages = [r["language"] for r in results]
    halley_iters = [r["halley_iter_mean"] for r in results]
    newton_iters = [r["newton_iter_mean"] for r in results]
    halley_us    = [r["halley_us"]        for r in results]
    newton_us    = [r["newton_us"]        for r in results]

    x = np.arange(len(languages))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    bars1 = ax.bar(x - width/2, halley_iters, width, label="Halley", color="#2E86AB")
    bars2 = ax.bar(x + width/2, newton_iters, width, label="Newton", color="#A23B72")
    ax.set_ylabel("Average Iterations", fontsize=12)
    ax.set_title(f"ProRail Benchmark: Average Iterations to Convergence\n"
                 f"(Halley vs Newton, N={bench['corpus_size']} real transitions)",
                 fontsize=13, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(languages, fontsize=10)
    ax.legend()
    ax.set_ylim(0, max(max(halley_iters), max(newton_iters)) * 1.4)
    ax.grid(axis="y", alpha=0.3)
    _annotate(ax, bars1, "{:.2f}")
    _annotate(ax, bars2, "{:.2f}")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "iterations_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()

    fig, ax = plt.subplots(figsize=(10, 6))
    bars1 = ax.bar(x - width/2, halley_us, width, label="Halley", color="#2E86AB")
    bars2 = ax.bar(x + width/2, newton_us, width, label="Newton", color="#A23B72")
    ax.set_ylabel("Median Time per Solve (μs)", fontsize=12)
    ax.set_title(f"ProRail Benchmark: Median Execution Time\n"
                 f"(Halley vs Newton across Languages, N={bench['corpus_size']})",
                 fontsize=13, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(languages, fontsize=10)
    ax.legend()
    ax.set_ylim(0, max(max(halley_us), max(newton_us)) * 1.4)
    ax.set_yscale("log")  # 50x dynamic range Python..C#
    ax.grid(axis="y", which="both", alpha=0.3)
    _annotate(ax, bars1, "{:.2f}")
    _annotate(ax, bars2, "{:.2f}")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "time_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()

    print("Regenerated iterations_comparison.png and time_comparison.png")


def _annotate(ax, bars, fmt: str) -> None:
    for bar in bars:
        h = bar.get_height()
        ax.annotate(fmt.format(h),
                    xy=(bar.get_x() + bar.get_width() / 2, h),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center", va="bottom", fontsize=9, fontweight="bold")


if __name__ == "__main__":
    main()
