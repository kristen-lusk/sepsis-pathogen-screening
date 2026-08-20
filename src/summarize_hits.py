#!/usr/bin/env python3
"""Collapse BLAST alignments to species-by-sample raw and normalized hit matrices."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def natural_key(text: str):
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", text)
    ]


def load_manifest(path: Path) -> tuple[list[str], dict[str, str]]:
    species_order = []
    accession_to_species = {}

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            accession = row["accession"].strip()
            species = row["species"].strip()
            accession_to_species[accession] = species
            species_order.append(species)

    return species_order, accession_to_species


def load_qc(path: Path) -> dict[str, int]:
    read_counts = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            read_counts[row["sample_id"]] = int(row["reads"])
    return read_counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "refs" / "reference_accessions.tsv",
    )
    parser.add_argument(
        "--qc",
        type=Path,
        default=ROOT / "results" / "sample_qc.tsv",
    )
    parser.add_argument(
        "--blast-dir",
        type=Path,
        default=ROOT / "results" / "blast",
    )
    parser.add_argument(
        "--raw-output",
        type=Path,
        default=ROOT / "results" / "pathogen_counts_by_sample.csv",
    )
    parser.add_argument(
        "--normalized-output",
        type=Path,
        default=ROOT / "results" / "pathogen_hits_per_100k_reads.csv",
    )
    args = parser.parse_args()

    species_order, accession_to_species = load_manifest(args.manifest)
    read_counts = load_qc(args.qc)

    blast_paths = sorted(
        args.blast_dir.glob("sample_*_vs_pathogen_panel.tsv"),
        key=lambda p: natural_key(p.name),
    )
    if not blast_paths:
        raise SystemExit(f"No BLAST TSV files found in {args.blast_dir}")

    counts = defaultdict(lambda: defaultdict(int))

    for path in blast_paths:
        sample_id = path.name.replace("_vs_pathogen_panel.tsv", "")
        seen_read_species = set()

        with path.open() as handle:
            for line_no, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue

                fields = line.split("\t")
                if len(fields) != 7:
                    raise SystemExit(
                        f"{path}:{line_no}: expected 7 BLAST columns, found {len(fields)}"
                    )

                qseqid, sseqid = fields[0], fields[1]
                if sseqid not in accession_to_species:
                    raise SystemExit(
                        f"{path}:{line_no}: subject accession {sseqid!r} "
                        "is not present in the reference manifest."
                    )

                species = accession_to_species[sseqid]
                key = (qseqid, species)
                if key in seen_read_species:
                    continue

                seen_read_species.add(key)
                counts[species][sample_id] += 1

    samples = sorted(read_counts, key=natural_key)
    missing_blast = [
        sample for sample in samples
        if not (args.blast_dir / f"{sample}_vs_pathogen_panel.tsv").exists()
    ]
    if missing_blast:
        raise SystemExit(
            "Missing BLAST result files for: " + ", ".join(missing_blast)
        )

    args.raw_output.parent.mkdir(parents=True, exist_ok=True)

    with args.raw_output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["species", *samples])
        for species in species_order:
            writer.writerow(
                [species, *[counts[species].get(sample, 0) for sample in samples]]
            )

    with args.normalized_output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["species", *samples])
        for species in species_order:
            row = [species]
            for sample in samples:
                total_reads = read_counts[sample]
                if total_reads <= 0:
                    raise SystemExit(f"Sample {sample} has zero reads in the QC table.")
                raw_count = counts[species].get(sample, 0)
                per_100k = (raw_count / total_reads) * 100000
                row.append(f"{per_100k:.3f}")
            writer.writerow(row)

    print(f"Wrote: {args.raw_output}")
    print(f"Wrote: {args.normalized_output}")


if __name__ == "__main__":
    main()
