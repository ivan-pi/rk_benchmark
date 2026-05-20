#!/usr/bin/env python3
"""
Generate a horizontal benchmark bar chart using matplotlib.

Usage (from repository root):
  python3 scripts/plot_mean_time_per_step.py
"""

from __future__ import annotations

import csv
import io
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT = REPO_ROOT / "scripts" / "mean_time_per_step_py.png"
AWK_SCRIPT = REPO_ROOT / "scripts" / "extract_mean_time_per_step.awk"


def extract_rows() -> tuple[list[str], list[float]]:
    with subprocess.Popen(
        ["./build/rk_benchmark"],
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        text=True,
    ) as bench_proc:
        result = subprocess.run(
            ["awk", "-f", str(AWK_SCRIPT)],
            stdin=bench_proc.stdout,
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        if bench_proc.stdout is not None:
            bench_proc.stdout.close()
        bench_proc.wait()
    if bench_proc.returncode != 0:
        raise RuntimeError("rk_benchmark exited with a non-zero status.")
    labels: list[str] = []
    values: list[float] = []
    reader = csv.reader(io.StringIO(result.stdout), delimiter="\t")
    for row in reader:
        if len(row) != 2:
            continue
        labels.append(row[0])
        values.append(float(row[1]))
    if not labels:
        raise RuntimeError("No benchmark rows parsed from rk_benchmark output.")
    return labels, values


def main() -> None:
    labels, values = extract_rows()
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
