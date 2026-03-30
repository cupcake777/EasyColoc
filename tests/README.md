# Tests

EasyColoc keeps tests lightweight and operational. The suite is designed to
answer three questions quickly:

1. Do the entry scripts still parse?
2. Do the main user-facing workflows still run on tiny synthetic inputs?
3. Do monitoring, reporting, bootstrap, and CLI helper paths still behave as expected?

## What To Run

| Goal | Command | Scope |
| --- | --- | --- |
| Fast syntax check | `Rscript tests/check_parse.R` | Parse-only validation for core R scripts |
| Standard local validation | `./easycoloc smoke` | End-to-end smoke suite used in routine repo checks |
| Single feature check | `Rscript tests/<file>.R` or `bash tests/<file>.sh` | Focused debugging for one subsystem |

## Test Map

| File | Covers | Included in `./easycoloc smoke` |
| --- | --- | --- |
| `check_parse.R` | Parse validation for core scripts | Yes |
| `smoke_test_plotting.R` | Locus plotting and synthetic panel rendering | Yes |
| `smoke_test_summary_report.R` | Result merging and HTML report generation | Yes |
| `smoke_test_runtime_monitor.R` | Runtime tracker, heartbeat, and monitor snapshot logic | Yes |
| `smoke_test_output_dir_tools.R` | Output-dir status, manifest, and monitoring helpers | Yes |
| `smoke_test_doctor.R` | Config and dependency doctor entrypoint | Yes |
| `smoke_test_refs.R` | Reference requirement listing | Yes |
| `smoke_test_bootstrap_refs.R` | Shared-reference bootstrap helpers | Yes |
| `smoke_test_bootstrap_demo.R` | Demo project bootstrap path | Yes |
| `smoke_test_1kg_setup.R` | Toy 1000 Genomes setup path | Yes |
| `smoke_test_gtex_bootstrap.R` | GTEx metadata bootstrap path | Yes |
| `smoke_test_ld_cache.R` | LD caching and reuse behavior | Yes |
| `smoke_test_cli.sh` | Shell entrypoints and `easycoloc init` scaffold | Yes |
| `smoke_test_plink_clump.R` | PLINK clumping behavior | No |

## Conventions

- Prefer synthetic or tiny fixture inputs over external downloads.
- Keep fixture build/population/version assumptions explicit in the fixture
  itself or in `tests/fixtures/README.md`.
- Keep each smoke test focused on one subsystem.
- Print a short `[SMOKE] ... passed` marker so failures are easy to localize in CI or terminal logs.
- Avoid writing outside temporary directories unless a test is explicitly exercising project-local behavior.
