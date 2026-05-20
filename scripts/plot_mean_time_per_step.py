#!/usr/bin/env python3
"""
Generate a horizontal benchmark bar chart using matplotlib.

Usage (from repository root):
  ./build/rk_benchmark | awk -f scripts/extract_mean_time_per_step.awk > results.txt
  python3 scripts/plot_mean_time_per_step.py results.txt
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT = REPO_ROOT / "scripts" / "mean_time_per_step_py.png"


def extract_rows(results_file: Path) -> tuple[list[str], list[float]]:
    labels: list[str] = []
    values: list[float] = []
    with results_file.open(newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) != 2:
                continue
            labels.append(row[0])
            values.append(float(row[1]))
    if not labels:
        raise RuntimeError(f"No benchmark rows parsed from: {results_file}")
    return labels, values


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python3 scripts/plot_mean_time_per_step.py <results.tsv>")

    results_file = Path(sys.argv[1]).expanduser()
    if not results_file.is_absolute():
        results_file = REPO_ROOT / results_file

    labels, values = extract_rows(results_file)
    fig, ax = plt.subplots(figsize=(12, 6.5))
    bars = ax.barh(labels, values, color="#6fa8dc")
    ax.set_xlabel("Mean time per step (us)")
    ax.set_ylabel("Approach")
    ax.set_title("RK23 Benchmark - Robertson Rate Equations")
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    ax.bar_label(bars, fmt="%.1f", padding=4)
    fig.tight_layout()
    fig.savefig(OUTPUT, dpi=150)
    print(f"Wrote plot: {OUTPUT}")


if __name__ == "__main__":
    main()
