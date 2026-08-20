#!/usr/bin/env python3
"""Download the reference accessions from NCBI and create one multi-FASTA panel."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def load_accessions(manifest: Path) -> list[str]:
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        accessions = [row["accession"].strip() for row in reader if row["accession"].strip()]

    if not accessions:
        raise SystemExit(f"No accessions found in {manifest}")
    return accessions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "refs" / "reference_accessions.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "refs" / "sepsis_pathogen_panel.fasta",
    )
    parser.add_argument(
        "--email",
        required=True,
        help="Email supplied to NCBI E-utilities as requested by NCBI guidance.",
    )
    args = parser.parse_args()

    accessions = load_accessions(args.manifest)

    params = {
        "db": "nuccore",
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
        "tool": "sepsis_pathogen_screening_portfolio",
        "email": args.email,
    }
    request = Request(
        f"{EFETCH_URL}?{urlencode(params)}",
        headers={"User-Agent": "sepsis-pathogen-screening-portfolio/1.0"},
    )

    print(f"Downloading {len(accessions)} NCBI nucleotide records...")
    with urlopen(request, timeout=120) as response:
        fasta_text = response.read().decode("utf-8")

    if not fasta_text.lstrip().startswith(">"):
        raise SystemExit("NCBI response did not appear to be FASTA data.")

    missing = [acc for acc in accessions if acc not in fasta_text]
    if missing:
        raise SystemExit(
            "Downloaded FASTA is missing expected accessions: " + ", ".join(missing)
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(fasta_text)
    print(f"Wrote reference panel: {args.output}")


if __name__ == "__main__":
    main()
