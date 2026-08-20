#!/usr/bin/env python3
"""Summarize FASTA record counts and sequence lengths for each input sample."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def natural_key(path: Path):
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", path.name)
    ]


def normalize_sample_id(stem: str) -> str:
    match = re.fullmatch(r"person_(\d+)", stem)
    return f"sample_{int(match.group(1))}" if match else stem


def qc_fasta(path: Path) -> tuple[int, int, float]:
    reads = 0
    total_bases = 0
    current_length = 0

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_length:
                    total_bases += current_length
                    current_length = 0
                reads += 1
            else:
                current_length += len(line)

    total_bases += current_length
    mean_length = total_bases / reads if reads else 0.0
    return reads, total_bases, mean_length


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT / "data" / "raw",
        help="Directory containing input FASTA files.",
    )
    parser.add_argument(
        "--pattern",
        default="person_*.fasta",
        help="Glob pattern for input FASTA files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "sample_qc.tsv",
    )
    args = parser.parse_args()

    fasta_paths = sorted(args.data_dir.glob(args.pattern), key=natural_key)
    if not fasta_paths:
        raise SystemExit(
            f"No FASTA files matched {args.pattern!r} in {args.data_dir}"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "reads", "total_bases", "mean_read_length"])

        for path in fasta_paths:
            reads, total_bases, mean_length = qc_fasta(path)
            writer.writerow(
                [
                    normalize_sample_id(path.stem),
                    reads,
                    total_bases,
                    f"{mean_length:.3f}",
                ]
            )

    print(f"Wrote QC summary: {args.output}")


if __name__ == "__main__":
    main()
