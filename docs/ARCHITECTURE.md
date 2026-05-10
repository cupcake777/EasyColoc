# Architecture Notes

This page is only for myself ^^

## Pipeline Flow

```text
GWAS input
  -> native harmonization
  -> PLINK locus discovery
  -> tabix QTL regional queries
  -> rsID/hash/position+allele matching
  -> ABF coloc
  -> optional SuSiE follow-up
  -> tables, plots, RDS bundles, runtime state, and reports
```

## Runtime Boundaries

- `easycoloc` is the user-facing CLI dispatcher.
- `run_coloc.R` runs the main pipeline.
- `src/` contains runtime R helpers loaded by the pipeline and tools.
- `tools/` contains CLI helpers for doctor, refs, bootstrap, reporting, QC, and monitoring.
- `tools/checks/` contains maintainer validation checks, not runtime pipeline code.
- `examples/checks/` contains user-facing example checks.

## Key Design Rules

- Keep genome build explicit for GWAS, QTL, LD, recombination, and annotation inputs.
- Keep large reference assets outside the repository.
- Prefer managed runs (`./easycoloc run --managed`) for long analyses because they write logs and runtime state consistently.
- Keep report-web generation as a read-only layer over completed result directories except for its own `report_web/` payload.

## Out Of Scope

- Inferring genome builds from ambiguous input files.
- Storing production reference panels in git.
- Treating smoke fixtures as biological evidence.
