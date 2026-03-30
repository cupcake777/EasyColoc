# EasyColoc

EasyColoc is a practical GWAS-to-QTL colocalization pipeline built for real
inputs rather than idealized examples. It focuses on two recurring pain points:

1. Public GWAS summary statistics often arrive in inconsistent `hg19`-era formats.
2. Many coloc workflows stop at `coloc.abf` and never carry signal-level follow-up through `SuSiE`.

EasyColoc standardizes GWAS inputs, queries tabix-indexed QTL resources, and
produces ABF and SuSiE results together with plots, manifests, runtime
snapshots, and an HTML report.

## At A Glance

| If you want to... | Run this |
| --- | --- |
| confirm the repo works locally | `./easycoloc smoke` |
| create the recommended local environment | `micromamba create -f environment.yml` |
| inspect required references and tools | `./easycoloc refs` and `./easycoloc doctor` |
| create a self-contained demo | `./easycoloc bootstrap-refs --demo ./demo_quickstart --run` |
| run the full pipeline with managed logging | `./easycoloc run --managed` |
| monitor an existing results directory | `./easycoloc status`, `./easycoloc monitor`, `./easycoloc watch` |

## Why EasyColoc

- `hg19` and `hg38` GWAS support with harmonization and fallback handling
- ABF + SuSiE in one pipeline instead of ABF-only coloc
- Build-aware 1000 Genomes bootstrap for `hg19` and `hg38`
- GTEx metadata bootstrap that can generate QTL summary CSVs and a ready-to-use YAML
- Managed runs with runtime heartbeat, manifest, monitor snapshot, and completion checks
- Demo mode that creates a tiny self-contained project and runs to HTML report

## 2-Minute Quickstart

### 1. Validate the installation

```bash
micromamba create -f environment.yml
micromamba activate easycoloc
./easycoloc refs
./easycoloc doctor
./easycoloc smoke
```

### 2. Create a self-contained demo

```bash
./easycoloc bootstrap-refs --demo ./demo_quickstart --run
```

This creates a tiny chr22 project, runs the pipeline end to end, and writes:

- `results/coloc_report.html`
- `results/all_colocalization_results.csv`
- `results/all_susie_results.csv`

### 3. Build a real reference panel

```bash
./easycoloc bootstrap-refs --setup-1kg ./refs/1kg_phase3_hg19 --build hg19 --pop EAS --chromosomes 1-22
./easycoloc bootstrap-refs --setup-1kg ./refs/1kg_phase3_hg38 --build hg38 --pop EAS --chromosomes 1-22
```

### 4. Prepare GTEx metadata

```bash
./easycoloc bootstrap-refs \
  --fetch-gtex-meta ./refs/gtex_meta \
  --gtex-eqtl-dir /path/to/gtex/eqtl \
  --gtex-sqtl-dir /path/to/gtex/sqtl
```

This downloads the GTEx sample attributes file, builds summary CSVs, and
generates `qtl_gtex_generated.yaml`.

## Typical Workflow

```text
validate environment
  -> bootstrap or point to references
  -> configure GWAS and QTL metadata
  -> run coloc
  -> inspect progress
  -> review merged results and report
```

## Architecture

```mermaid
flowchart TD
    A[GWAS summary statistics<br/>hg19 or hg38] --> B[Harmonization layer<br/>gwaslab plus fallback]
    B --> C[Locus discovery<br/>PLINK clump]
    C --> D[QTL sigPairs prefilter<br/>skip empty datasets fast]
    C --> E[QTL allPairs regional query]
    D --> F[Variant matching<br/>rsID or hash or allele]
    E --> F
    F --> G[Coloc engine<br/>ABF plus SuSiE]
    G --> H[CSV summaries]
    G --> I[Locus plots]
    G --> J[RDS bundles]
    G --> K[HTML report]
```

See [ARCHITECTURE.md](docs/ARCHITECTURE.md) for a fuller description.

## Main Commands

| Command | Use it for |
| --- | --- |
| `./easycoloc run --managed` | Full pipeline run with managed logging and completion checks |
| `./easycoloc check /path/to/output_dir` | Determine whether a run finished cleanly |
| `./easycoloc status /path/to/output_dir` | Summarize per-GWAS progress and output counts |
| `./easycoloc monitor /path/to/output_dir` | Print the latest runtime heartbeat and output snapshot |
| `./easycoloc watch /path/to/output_dir 60 logs/monitor/my_run.log` | Append periodic monitor snapshots to a local log |
| `./easycoloc manifest /path/to/output_dir` | Build an output manifest for an existing results directory |
| `./easycoloc refs --include-qtl-files` | Inspect required and optional references |
| `./easycoloc bootstrap-refs ...` | Materialize local references or demo assets |
| `./easycoloc smoke` | Run the standard local smoke suite |

## Repository Layout

- `src/`: core R modules
- `tools/`: executable helper scripts
- `tests/`: smoke and regression tests
- `config/`: portable public defaults plus generated examples
- `data/`: optional local input staging area for ad hoc runs
- `docs/`: user-facing documentation
- `examples/`: minimal demos
- `templates/`: `easycoloc init` scaffold

Detailed layout notes are in [REPO_LAYOUT.md](docs/REPO_LAYOUT.md).

## Key Outputs

| Output | Meaning |
| --- | --- |
| `all_colocalization_results.csv` | merged locus-level ABF coloc results |
| `significant_colocalizations_PP4_*.csv` | thresholded coloc hits |
| `all_susie_results.csv` | merged SuSiE output |
| `plots/*.pdf|png` | locus plots |
| `rds/*.rds` | serialized locus bundles for rerendering or inspection |
| `coloc_report.html` | interactive report |
| `output_manifest.tsv` | machine-readable output inventory |

## Documentation

- [TUTORIAL.md](docs/TUTORIAL.md)
- [ARCHITECTURE.md](docs/ARCHITECTURE.md)
- [REFERENCE_DATA.md](docs/REFERENCE_DATA.md)
- [DOCKER.md](docs/DOCKER.md)
- [COMPETITIVE_POSITIONING.md](docs/COMPETITIVE_POSITIONING.md)

## Validation

Run the standard local validation suite:

```bash
./easycoloc smoke
```

For a lighter parse-only check:

```bash
Rscript tests/check_parse.R
```

For a test inventory and per-file coverage notes, see
[tests/README.md](tests/README.md).

For config layering and local-private override guidance, see
[config/README.md](config/README.md).

To monitor an existing results directory without writing into that directory:

```bash
./easycoloc watch /path/to/output_dir 60
```
