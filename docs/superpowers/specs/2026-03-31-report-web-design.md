# EasyColoc Web Report And Docs Restructure Design

Date: 2026-03-31
Status: Proposed
Scope: design only, no implementation in this document

## 1. Problem Statement

EasyColoc currently generates an HTML report from `run_coloc.r` through
`src/utils_report.R`, but the report layer is tightly coupled to string-built R
HTML and is not well positioned for richer interaction, clearer navigation, or
maintainable visual polish.

The documentation has a second problem: the current entry points are split
across `README.md`, `docs/TUTORIAL.md`, and several supporting docs. This makes
the "what should I do first?" path slower than necessary, especially for users
who already have a completed `results/` directory and mainly want to inspect it.

## 2. User Goal

The primary target user for this change is:

- a user who already has a completed EasyColoc `results/` directory
- wants to inspect that run locally
- wants one CLI command to launch a polished interactive web report

Secondary support remains valuable for:

- new users validating the repo with a demo
- users who need architecture or reference documentation

## 3. Goals

- Add a local web-based interactive report for an existing `results/` directory
- Keep `./easycoloc` as the only required user-facing entry point
- Allow lightweight data transformation at report startup so the frontend can
  consume stable JSON instead of raw heterogeneous CSV directly
- Restructure documentation so the shortest path for existing-result users is
  obvious and authoritative
- Remove or demote docs that do not materially help users start, inspect
  results, or understand the system

## 4. Non-Goals

- No hosted SaaS or multi-user report service
- No requirement to support static-file sharing as the primary delivery mode
- No rewrite of the analysis pipeline itself
- No attempt to make the frontend the source of truth for scientific semantics
- No speculative biological interpretation beyond existing run outputs

## 5. Chosen Approach

Selected approach: `frontend sub-application + report data adapter + CLI
orchestration`

This corresponds to:

- keeping EasyColoc analysis outputs as the source of truth
- introducing a dedicated `web/` frontend application for interaction and
  visual presentation
- introducing a lightweight transformation step that converts `results/`
  artifacts into a stable JSON payload for the frontend
- exposing the whole experience through a new CLI command

This approach is preferred over extending the current single-file HTML report
because:

- R string-built HTML will become harder to maintain as interaction grows
- frontend styling, componentization, and navigation are cleaner in a dedicated
  web stack
- a thin data adapter keeps scientific output semantics anchored in the current
  EasyColoc pipeline instead of duplicating analysis logic in JavaScript

## 6. Proposed Architecture

The report system should be split into three layers.

### 6.1 Results Layer

Existing pipeline outputs remain unchanged in role:

- `all_colocalization_results.csv`
- `all_susie_results.csv`
- `significant_colocalizations_PP4_*.csv`
- plot assets
- runtime tracker files
- manifest and summary artifacts

These files remain the authoritative outputs of the analysis run.

### 6.2 Report Data Adapter Layer

A lightweight adapter generates frontend-ready JSON from a chosen `results/`
directory at report launch time.

Responsibilities:

- validate required files and detect partial runs
- load summary tables and runtime artifacts
- normalize field names needed by the frontend
- derive compact summary metrics
- index plots and downloadable assets
- write a stable JSON bundle for the frontend

Constraints:

- must not invent results not present in current outputs
- must explicitly encode missing or failed artifacts
- should remain lightweight enough to run on report launch

### 6.3 Web Frontend Layer

A dedicated frontend application under `web/` renders the report UI from the
adapter JSON.

Responsibilities:

- overview dashboard
- interactive filtering and sorting
- drill-down into a single locus result
- asset browsing and audit-oriented run inspection

The frontend should not parse arbitrary EasyColoc raw outputs directly when a
stable adapted payload exists.

### 6.4 CLI Orchestration Layer

Add a new command:

```bash
./easycoloc report-web /path/to/results
```

The CLI will:

1. validate the provided results directory
2. run or refresh the report data adapter
3. start the local web service
4. print the access URL
5. optionally open the browser if requested

## 7. CLI Design

### 7.1 Command

```bash
./easycoloc report-web /path/to/results
```

### 7.2 Proposed Flags

- `--host`
- `--port`
- `--open`
- `--no-open`
- `--refresh-data`
- `--project-name`

### 7.3 Behavioral Contract

- If the directory is missing, fail fast with a clear path-specific error
- If key outputs are missing, start only if the available data are still
  sufficient for a partial report, and label the report as partial
- If transformed report data already exist, reuse them unless
  `--refresh-data` is supplied
- The command must not require the user to manually enter the `web/` directory
  or run `npm` commands directly

## 8. Frontend Information Architecture

The report is for inspecting one completed run, not for generic marketing or
project landing content.

### 8.1 Overview

Purpose: answer "is this run worth deeper inspection?" in seconds.

Should include:

- project name
- result directory
- report generation time
- total tests
- PP4 summary metrics
- count of strong or suggestive hits
- SuSiE availability summary
- run integrity or warning badges
- preview of top hits

### 8.2 Explore Results

This is the primary working surface.

Should support:

- search across key identifiers
- filtering by trait, QTL dataset, tissue, gene or phenotype, and PP4 ranges
- sortable result table
- selectable rows with linked detail view
- visible handling of missing fields

This page should replace the current workflow of manually opening CSV files and
hunting for related plots or summary values.

### 8.3 Locus Detail

Purpose: inspect a single result without leaving the report context.

Should include:

- locus metadata
- coloc metrics
- SuSiE summary when available
- source file provenance
- links or previews for plots and related assets
- clear markers for missing outputs or skipped steps

### 8.4 Run Assets

Purpose: support handoff, auditing, and reproducibility review.

Should include:

- output inventory
- manifest summary
- runtime or heartbeat snapshot summary
- key generated files
- download links where appropriate

## 9. Data Contract For The Frontend

The frontend should consume a stable report payload instead of raw output files.

Suggested logical sections:

- `meta`
- `summary`
- `results`
- `susie_summary`
- `assets`
- `runtime`
- `warnings`

Notes:

- `meta` holds project name, result path, generation time, and report version
- `summary` holds overview counters and headline metrics
- `results` holds normalized row-level coloc entries used by the main explorer
- `susie_summary` holds optional fine-mapping-specific compact records
- `assets` holds indexed plot and downloadable file locations
- `runtime` holds monitor or completion state relevant for auditing
- `warnings` holds explicit partial-run or missing-artifact flags

## 10. Visual And Product Direction

The report should look like a serious scientific results browser, not a generic
admin dashboard.

Design direction:

- clean but distinctive typography and spacing
- strong visual hierarchy around summary and result exploration
- restrained but intentional color system
- explicit warning or incomplete-state treatments
- responsive layout for laptop-first usage, with workable mobile fallback

The first screen should prioritize interpretation efficiency over decorative
cards.

## 11. Documentation Restructure

### 11.1 README.md

`README.md` should become the authoritative top-level entry point.

Its structure should prioritize:

1. what EasyColoc is
2. the fastest path for users with an existing results directory
3. the minimal demo path for first-time validation
4. links to deeper docs for configuration and architecture

The README should stop trying to be the complete tutorial and should instead
route users by intent.

### 11.2 docs/TUTORIAL.md

`docs/TUTORIAL.md` should be reorganized by user scenario:

- existing results directory -> open the interactive report
- demo run -> generate a tiny example and inspect it
- full real analysis -> configure references and run the pipeline

This is better than a single long narrative because the primary user journeys
are different.

### 11.3 docs/ARCHITECTURE.md

`docs/ARCHITECTURE.md` should be kept, but rewritten more directly where
necessary so it explains implementation boundaries rather than project rhetoric.

It should gain a short section describing:

- report data adapter
- web frontend
- CLI orchestration for report launch

### 11.4 Docs To Demote Or Remove

`docs/COMPETITIVE_POSITIONING.md` should be reviewed first for removal or
demotion because it does not directly support startup, results inspection, or
system understanding for most users.

General removal rule:

- remove or demote docs that duplicate command instructions already covered by
  `README.md` and `docs/TUTORIAL.md`
- remove or demote docs that contain internal positioning language without
  helping a user complete a task

## 12. Fast-Start Experience

The primary fast-start path after this redesign should be:

```bash
./easycoloc report-web /path/to/results
```

The user should be able to understand what happens next from the terminal
output alone:

- whether the results directory was accepted
- whether the adapter data were generated or refreshed
- which local URL to open
- whether the run appears complete or partial

The documentation should mirror this exact path and not require users to infer
it from longer pipeline material.

## 13. Risks And Mitigations

### Risk 1: frontend and analysis semantics drift apart

Mitigation:

- keep the frontend on a stable adapted payload
- keep scientific field interpretation in the adapter layer
- do not reimplement analysis semantics in client components

### Risk 2: report launch becomes too heavy

Mitigation:

- keep transformation lightweight
- cache adapted report data where possible
- refresh only on explicit request or stale detection

### Risk 3: documentation becomes fragmented again

Mitigation:

- establish `README.md` as the single top-level entry point
- reduce repeated command descriptions in secondary docs
- route by scenario instead of duplicating full command sequences everywhere

### Risk 4: partial or failed runs are presented as healthy

Mitigation:

- surface missing artifacts and runtime failures explicitly
- add visible warnings in both CLI output and the web UI
- never silently substitute absent outputs with normal-state placeholders

## 14. Implementation Boundaries For The Next Phase

This design implies the following implementation work, but does not perform it:

- add a frontend app under `web/`
- add a report data adapter command or module
- add `report-web` to the CLI
- revise report-related docs
- audit candidate docs for deletion or demotion
- add verification for the new report-launch path

## 15. Acceptance Criteria

The redesign is successful when all of the following are true:

- a user with an existing results directory can launch the report with one CLI
  command
- the report provides overview, results exploration, locus detail, and asset
  audit views
- the frontend reads stable adapted report data rather than arbitrary raw
  outputs directly
- the README makes the existing-results fast path obvious
- tutorial and architecture docs no longer duplicate or obscure the primary
  usage path
- non-essential or low-value docs are removed or clearly demoted

