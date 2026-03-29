# Competitive Positioning

## Goal

EasyColoc should compete on end-to-end coloc execution quality, not on marketing
claims. The product story is strongest when it does three things better than
the current alternatives:

1. Launch faster than `ColocQuiaL`
2. Operationalize better than an ad hoc `gwaslab + custom scripts` stack
3. Produce outputs that are easier to audit, resume, and hand off

## Honest Comparison

| Capability | gwaslab | ColocQuiaL | EasyColoc |
|---|---|---|---|
| GWAS munging / harmonization | Excellent | Limited | Good via gwaslab integration |
| QTL-GWAS batch coloc | Partial | Strong | Strong |
| Unified CLI | Strong | Weak | Strong |
| Portable config paths | Strong | Weak | Strong |
| Managed long runs | Limited | Limited | Strong |
| Runtime heartbeat / task checkpointing | Limited | Limited | Strong |
| Built-in output manifest / run status | Limited | Limited | Strong |
| Minimal project scaffold | Moderate | Weak | Strong |
| Publication-ready coloc plotting | Moderate | Moderate | Strong |

## What EasyColoc Should Emphasize

- `gwaslab` is a dependency and force multiplier, not an enemy. EasyColoc wins by making `gwaslab` harmonization the beginning of the workflow rather than the end.
- `ColocQuiaL` is the more direct coloc competitor. EasyColoc should consistently beat it on setup friction, monitoring, resume safety, and reportability.
- The repo should market itself as a reproducible coloc platform rather than just a script collection.

## High-Impact Directions

### 1. Product Surface

- Keep `./easycoloc` as the canonical entrypoint.
- Make `doctor`, `init`, `status`, `monitor`, and `manifest` first-class and stable.
- Prefer portable config-relative paths over hard-coded repository assumptions.

### 2. Proof, Not Claims

- Add benchmark datasets with fixed expected output counts and wall-clock summaries.
- Publish a reproducible comparison notebook against one or two public loci.
- Keep smoke tests and lightweight CI green at all times.

### 3. Star Drivers

- A landing page that demonstrates one-command setup and one-command health checks
- A tutorial that goes from public GWAS summary statistics to coloc report in under 10 minutes
- Clear examples for eQTL, sQTL, pQTL, and APA-style phenotype identifiers
- Honest scope boundaries: what EasyColoc does, what it delegates to gwaslab, and what remains user-supplied

## Near-Term Execution Plan

1. Keep reducing hard-coded paths and repository-root assumptions.
2. Add benchmarkable public examples with expected output manifests.
3. Add lightweight CI for parse checks, shell syntax, and doctor/template validation.
4. Add API-like machine-readable summaries so dashboards and notebooks can inspect runs without scraping logs.
