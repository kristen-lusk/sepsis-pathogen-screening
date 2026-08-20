# Input data

The original course project used ten blood-derived sequencing FASTA files named
`person_1.fasta` through `person_10.fasta`, with 100,000 reads per sample.

Those raw sequencing files are not distributed in this public portfolio repository.
Redistribution permission/provenance for the course-provided sequence data was not
established during portfolio preparation, so the repository intentionally contains only
code, reference accessions, compact summary results, and figures.

To reproduce the analysis with authorized copies of the original inputs, place them in:

```text
data/raw/
├── person_1.fasta
├── person_2.fasta
├── ...
└── person_10.fasta
```

The public-facing scripts convert the original `person_N` filenames to `sample_N`
identifiers in generated tables and figures.

Do not commit raw sequencing files to the repository unless you have confirmed that
you are permitted to redistribute them.
