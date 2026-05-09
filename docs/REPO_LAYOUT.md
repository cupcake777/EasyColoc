# Repository Layout

This repository separates source code, user configuration, documentation,
tests, and local run artifacts. The root should stay small: one CLI, one main
pipeline entrypoint, and project metadata.

## Reading Order

1. `README.md` for the short product and usage overview
2. `docs/TUTORIAL.md` for setup and execution
3. `config/README.md` for public vs private config layers
4. `docs/ARCHITECTURE.md` for pipeline internals
5. `tests/README.md` for validation coverage

## Root Files

| Path | Purpose |
| --- | --- |
| `easycoloc` | Unified shell CLI |
| `run_coloc.r` | Main R pipeline entrypoint |
| `README.md` | Project overview and common commands |
| `environment.yml` | Recommended micromamba environment |
| `Dockerfile` | Container build recipe |
| `CITATION.cff` | Citation metadata |
| `LICENSE`, `CONTRIBUTING.md` | Project metadata |

## Source And Operations

| Path | Purpose |
| --- | --- |
| `src/` | R modules for config, formatting, harmonization, coloc, plotting, reporting, references, bootstrap, and runtime state |
| `tools/` | Operational scripts called by the CLI or maintainers |
| `web/` | Report web UI source; `web/dist/` and `web/node_modules/` are local build/install artifacts |

Important tool groups:

| Path | Purpose |
| --- | --- |
| `tools/doctor_easycoloc.R` | Validate configs, required files, and external tools |
| `tools/list_reference_requirements.R` | Print reference inventory |
| `tools/bootstrap_references.R` | Create demo projects, bootstrap 1KG, or prepare GTEx metadata |
| `tools/build_harmony_qc_report.R` | QC harmonized GWAS cache files |
| `tools/run_pipeline_managed.sh` | Run pipeline with managed logs and post-run checks |
| `tools/monitor_easycoloc.R`, `tools/summarize_run_status.R`, `tools/check_run_completion.R` | Inspect existing result directories |
| `tools/build_report_web_data.R`, `tools/start_report_web.sh`, `tools/report_web_server.mjs` | Build and serve the local report UI |

## Configuration And Templates

| Path | Purpose |
| --- | --- |
| `config/global.yaml` | Portable default global settings |
| `config/gwas.yaml` | Example GWAS dataset config |
| `config/qtl.yaml` | Example QTL config |
| `config/QTL_summary.csv` | Example QTL metadata table used by `config/qtl.yaml` |
| `config/qtl_gtex_generated.yaml` | Example generated GTEx-style QTL config |
| `config/local/README.md` | Guidance for private local configs |
| `templates/project/` | Files copied by `./easycoloc init TARGET_DIR` |

`config/local/*` is ignored except for its README. Put lab-specific absolute
paths there instead of editing public example configs.

## Tests And Examples

| Path | Purpose |
| --- | --- |
| `tests/check_parse.R` | Parse all important R scripts |
| `tests/smoke_test_*.R`, `tests/smoke_test_*.sh` | Focused smoke/regression checks |
| `tests/fixtures/` | Small coordinate-compatible fixtures |
| `tests/manual/` | Manual diagnostics outside the standard smoke suite |
| `examples/minimal/` | Minimal synthetic plotting demo |

Run the standard validation suite with:

```bash
./easycoloc smoke
```

## Documentation

| Path | Purpose |
| --- | --- |
| `docs/TUTORIAL.md` | End-to-end usage tutorial |
| `docs/ARCHITECTURE.md` | Pipeline architecture |
| `docs/REFERENCE_DATA.md` | Reference setup and bootstrap notes |
| `docs/REFERENCE_COMPATIBILITY.md` | Build/population compatibility guidance |
| `docs/DOCKER.md` | Container usage |
| `docs/assets/` | README/tutorial assets |

`docs/notes/`, `docs/optimization/`, and `docs/superpowers/` are treated as
local development notes by `.gitignore`. They may exist in a working directory
but are not part of the public operational surface.

## Local Run Artifacts

The following paths are ignored and can be deleted/recreated locally:

| Path | Contents |
| --- | --- |
| `results/` | Pipeline outputs, reports, manifests, runtime state |
| `harmony/` | Reusable harmonized GWAS caches |
| `temp/` | PLINK, tabix, and bootstrap intermediates |
| `logs/` | Managed-run logs and monitor logs |
| `examples/minimal/output/` | Generated demo figures |
| `web/dist/`, `web/node_modules/` | Web build and dependency artifacts |

To inspect the clean tracked structure only:

```bash
git ls-files
```

To find ignored local clutter:

```bash
git status --ignored --short
```
