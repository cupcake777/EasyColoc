# Architecture

## Design Goal

EasyColoc is meant to make real coloc work operational, not just statistically
possible. The pipeline is therefore organized around four stages:

1. Normalize messy GWAS inputs
2. Minimize unnecessary QTL scanning
3. Run both ABF and SuSiE coloc
4. Leave behind outputs that are easy to inspect and hand off

## Pipeline Overview

```mermaid
flowchart LR
    A[GWAS Input] --> B[Harmonization]
    B --> C[Locus Discovery]
    C --> D[QTL Regional Query]
    D --> E[Variant Matching]
    E --> F[ABF Coloc]
    F --> G[SuSiE Coloc]
    G --> H[Summaries]
    G --> I[Plots]
    G --> J[RDS Bundles]
    G --> K[HTML Report]
    H --> L[Report Web Payload]
    I --> L
    J --> L
    L --> M[Local Report Web Server]
```

## Main Components

### 1. Harmonization

- Implemented in `run_coloc.r` and `src/utils_format.R`
- Uses native EasyColoc harmonization when reference assets are available
- Supports fallback handling for liftOver-related failure paths
- Keeps build-aware reference selection explicit

### 2. Locus Discovery

- Uses PLINK clumping on a build-matched 1000 Genomes panel
- Build-aware PLINK reference selection now supports both `plink_hg19` and `plink_hg38`
- Population-specific keep files are reused for both clumping and SuSiE LD

### 3. QTL Query Layer

- QTL data are expected as tabix-indexed allPairs/sigPairs files
- `sigPairs` prefiltering is enabled by default to avoid expensive allPairs scans when a dataset has no local signal
- Relative QTL metadata paths are resolved robustly for portable projects

### 4. Variant Matching

EasyColoc uses a three-tier match strategy:

1. direct rsID overlap
2. rsID rescue through dbSNP hash tables
3. position + allele matching

This layer is the bridge between inconsistent GWAS inputs and consortium-style
QTL files.

### 5. Coloc Engine

- ABF coloc is always attempted
- SuSiE coloc is triggered when signal strength and LD prerequisites are met
- Credible set extraction and plotting hooks are integrated into the same run

### 6. Output Layer

| Output type | Purpose |
| --- | --- |
| merged CSV summaries | downstream analysis and filtering |
| SuSiE tables | fine-mapping review |
| locus plots | visual interpretation and handoff |
| serialized RDS bundles | rerendering and debugging |
| runtime tracker files | monitoring and auditing |
| HTML report | compact review artifact |
| output manifest | machine-readable output inventory |

## Bootstrap Layer

Bootstrap functionality is a first-class part of the architecture:

- `--demo`: creates a toy project and can run it end to end
- `--setup-1kg`: downloads build-specific 1000 Genomes VCFs and converts them to PLINK
- `--fetch-gtex-meta`: downloads GTEx sample attributes and can generate QTL summary metadata plus YAML

These paths live primarily in:

- `tools/bootstrap_references.R`
- `src/utils_bootstrap.R`
- `src/utils_download.R`

## Runtime And Observability

EasyColoc tracks long runs explicitly:

- `active_run.json`
- `heartbeat.json`
- `task_state.tsv`
- `event_log.ndjson`
- `monitor_snapshot.json`

This keeps the pipeline inspectable without scraping unstable console output.

## Report Web Layer

The report web path is a lightweight layer on top of finished outputs so users
can inspect existing result directories without rerunning colocalization.

### 1. Adapter

- CLI command: `./easycoloc report-web /path/to/results`
- Entrypoint scripts: `tools/build_report_web_data.R`, `tools/start_report_web.sh`
- Role: validate output directory shape and map EasyColoc outputs into web-ready inputs

### 2. Payload

- Generated artifact: `report_web/report-data.json` under the target results directory
- Location: inside the `report_web/` subdirectory of the selected results directory
- Role: normalized payload for the frontend, built from summary tables and report assets

### 3. Local Server

- Server script: `tools/report_web_server.mjs`
- Role: serve static web assets and payload locally for interactive browsing

### 4. CLI Wiring

- User-facing entrypoint: `easycoloc` dispatcher (`report-web` subcommand)
- Role: run payload generation first, then launch local server with the selected results directory

## What Is Deliberately Not In Scope

- full GWAS munging coverage for common summary-statistic formats
- storage of large reference assets inside the repository
- implicit magic around build inference when explicit configuration is safer
