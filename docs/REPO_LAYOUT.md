# Repository Layout

EasyColoc keeps user-facing entrypoints shallow and pushes implementation into
responsibility-based directories.

## Reading Order

If you are new to the repo, start here:

1. `README.md` for the product view
2. `docs/TUTORIAL.md` for setup and execution
3. `docs/ARCHITECTURE.md` for pipeline structure
4. `config/README.md` for config layering
5. `tests/README.md` for validation entrypoints

## Root

Only true entrypoints and project metadata should remain at the repository
root:

- `easycoloc`: unified CLI
- `run_coloc.r`: main pipeline entrypoint
- `README.md`: project overview
- `Dockerfile`, `LICENSE`, `CONTRIBUTING.md`

## Primary Directories

| Directory | Purpose | Typical reader |
| --- | --- | --- |
| `src/` | reusable R modules for config, bootstrap, harmonization, coloc, plotting, output, and runtime tracking | developers changing pipeline behavior |
| `tools/` | executable helpers such as bootstrap, doctor, monitor, watch, manifest, and plot reruns | users and maintainers |
| `tests/` | parse checks and smoke tests for core workflows | anyone validating a change |
| `config/` | portable defaults, local-private config guidance, and example/generated metadata tables | users adapting the pipeline |
| `data/` | optional local staging area for repo-local ad hoc inputs | users experimenting locally |
| `docs/` | tutorial, architecture, reference setup, and operational notes | users and reviewers |
| `examples/` | small demos and lightweight outputs | first-time users |
| `templates/` | files used by `easycoloc init` | users creating portable projects |

## Local-Only / Derived Content

The following locations are intentionally treated as local working state and are
ignored by git:

| Path | Why it stays local |
| --- | --- |
| `logs/` | managed-run logs, watch logs, and local debugging output |
| `harmony/` | cached harmonized GWAS tables |
| `temp/` | transient PLINK and bootstrap intermediates |
| `results/` | repository-local outputs from ad hoc runs |
| `docs/notes/`, `docs/optimization/` | optional local development notes |

## Current Placement Notes

| Path | Note |
| --- | --- |
| `config/README.md` | explains public defaults vs local-private overrides |
| `config/QTL_summary.csv` | default QTL metadata table consumed by `config/qtl.yaml` |
| `config/qtl_gtex_generated.yaml` | generated GTEx-oriented QTL config example |
| `tools/rerun_plots.R` | rerender locus plots from existing RDS bundles without rerunning coloc |
| `tools/watch_output_dir.sh` | append periodic monitor snapshots for an arbitrary results directory to a local log |
| `tests/manual/` | one-off diagnostics that are intentionally outside the standard smoke suite |
