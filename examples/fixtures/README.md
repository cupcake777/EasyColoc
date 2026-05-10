# Fixture Design Notes

The fixtures under `examples/fixtures/` are not arbitrary placeholders. They are
small, synthetic inputs chosen to preserve the compatibility constraints that
matter for scientific interpretation.

## Plotting Fixture

The plotting smoke test uses:

- coordinate system: `hg38`
- chromosome: `chr1`
- phenotype anchor: `ENST00000654683.1|CCDC30|chr1:42482663-42484158|+`
- recombination population: `CHB`

Files:

- `annotation/smoke_hg38_chr1.gtf`
- `recomb/hg38/CHB/CHB_recombination_map_hapmap_format_hg38_chr_1.txt`

## Why These Match

- The GTF fixture places `CCDC30` at the same `hg38 chr1` interval used by the
  synthetic phenotype label in `tools/checks/smoke_test_plotting.R`.
- The recombination fixture uses the repository's expected `hg38` HapMap-style
  naming convention and explicitly encodes the `CHB` population in the path and
  filename.
- The fixture is intentionally minimal, but it is still build-aware and
  population-explicit. That prevents the smoke figure from looking visually
  complete while being scientifically inconsistent.

## What This Fixture Is Not

- It is not a production annotation release.
- It is not intended for biological interpretation.
- It is not evidence that `CHB` is the correct population for a real user run.

Its purpose is narrower: make the plotting smoke test reproducible while
preserving the compatibility logic that a real analysis should also follow.
