#!/usr/bin/env python3
"""Generate the two portfolio figures from the summarized count matrix."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]


def natural_key(text: str):
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", text)
    ]


def load_counts(path: Path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    samples = sorted(
        [field for field in reader.fieldnames if field != "species"],
        key=natural_key,
    )
    return samples, rows


def main() -> None:
    input_path = ROOT / "results" / "pathogen_counts_by_sample.csv"
    output_dir = ROOT / "results" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    samples, rows = load_counts(input_path)

    candida_row = next(
        row for row in rows if row["species"] == "Candida albicans"
    )
    candida_counts = [int(candida_row[sample]) for sample in samples]

    plt.figure(figsize=(9, 5))
    plt.bar(samples, candida_counts)
    plt.ylabel("Unique read-species hits")
    plt.xlabel("Sample")
    plt.title("Candida albicans alignment signal by sample")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_dir / "candida_hits_by_sample.png", dpi=180)
    plt.close()

    species = [row["species"] for row in rows]
    matrix = [
        [math.log10(int(row[sample]) + 1) for sample in samples]
        for row in rows
    ]

    plt.figure(figsize=(10, 7))
    image = plt.imshow(matrix, aspect="auto")
    plt.colorbar(image, label="log10(unique hits + 1)")
    plt.xticks(range(len(samples)), samples, rotation=45, ha="right")
    plt.yticks(range(len(species)), species)
    plt.xlabel("Sample")
    plt.ylabel("Reference organism")
    plt.title("Targeted pathogen-screening alignment signals")
    plt.tight_layout()
    plt.savefig(output_dir / "pathogen_hits_heatmap.png", dpi=180)
    plt.close()

    print(f"Wrote figures to {output_dir}")


if __name__ == "__main__":
    main()
