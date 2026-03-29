# Repository Layout

EasyColoc keeps user-facing entrypoints shallow and pushes implementation into
responsibility-based directories.

## Root

Only true entrypoints and project metadata should remain at the repository
root:

- `easycoloc`: unified CLI
- `run_coloc.r`: main pipeline entrypoint
- `README.md`: project overview
- `Dockerfile`, `LICENSE`, `CONTRIBUTING.md`

## Primary Directories

- `src/`: reusable R modules for config, bootstrap, harmonization, coloc, plotting, output, and runtime tracking
- `tools/`: executable helper scripts such as bootstrap, doctor, monitor, manifest, plotting reruns, and conversion utilities
- `tests/`: smoke tests, parse checks, and manual diagnostic scripts
- `config/`: default pipeline configs and example/generated QTL metadata tables
- `docs/`: user documentation, architecture notes, tutorial, and reference setup guides
- `examples/`: minimal demos and lightweight example outputs
- `templates/`: files used by `easycoloc init`

## Local-Only / Derived Content

The following locations are intentionally treated as local working state and are
ignored by git:

- `logs/`: managed-run logs and local debugging output
- `harmony/`: cached harmonized GWAS tables
- `temp/`: transient PLINK and bootstrap intermediates
- `results/`: repository-local outputs from ad hoc runs
- `docs/notes/`, `docs/optimization/`: development notes kept locally when needed

## Current Placement Notes

- `config/QTL_summary.csv`: default QTL metadata table consumed by `config/qtl.yaml`
- `config/qtl_gtex_generated.yaml`: generated GTEx-oriented QTL config example
- `tools/rerun_plots.R`: rerender locus plots from existing RDS bundles without rerunning the full pipeline
- `tests/manual/`: one-off developer diagnostics that are useful locally but are not part of the standard smoke suite
