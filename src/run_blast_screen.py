#!/usr/bin/env python3
"""Run BLASTn for each sequencing sample against the targeted pathogen panel."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT / "data" / "raw",
    )
    parser.add_argument(
        "--pattern",
        default="person_*.fasta",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=ROOT / "refs" / "sepsis_pathogen_panel.fasta",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=ROOT / "results" / "blast",
    )
    args = parser.parse_args()

    if shutil.which("blastn") is None:
        raise SystemExit(
            "blastn was not found on PATH. Install NCBI BLAST+ before running this step."
        )

    if not args.reference.is_file() or args.reference.stat().st_size == 0:
        raise SystemExit(
            f"Reference FASTA is missing or empty: {args.reference}\n"
            "Run prepare_reference.py first."
        )

    fasta_paths = sorted(args.data_dir.glob(args.pattern), key=natural_key)
    if not fasta_paths:
        raise SystemExit(
            f"No FASTA files matched {args.pattern!r} in {args.data_dir}"
        )

    args.results_dir.mkdir(parents=True, exist_ok=True)

    for fasta_path in fasta_paths:
        sample_id = normalize_sample_id(fasta_path.stem)
        out_tsv = args.results_dir / f"{sample_id}_vs_pathogen_panel.tsv"

        cmd = [
            "blastn",
            "-query", str(fasta_path),
            "-subject", str(args.reference),
            "-task", "blastn",
            "-outfmt", "6 qseqid sseqid pident length qlen evalue bitscore",
            "-max_target_seqs", "5",
            "-evalue", "1e-5",
            "-word_size", "11",
            "-out", str(out_tsv),
        ]

        print(f"Running BLASTn for {sample_id}...")
        subprocess.run(cmd, check=True)

    print(f"BLAST outputs written to {args.results_dir}")


if __name__ == "__main__":
    main()
