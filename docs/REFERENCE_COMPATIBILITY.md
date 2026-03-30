# Reference Compatibility

This document makes the build, population, and annotation compatibility rules
explicit. For colocalization workflows, the most common source of silent
scientific error is not software failure but a reference mismatch that still
allows the pipeline to run.

## Core Principle

For every run, track these dimensions together:

- GWAS genome build
- QTL genome build
- LD reference panel build
- LD reference panel population
- recombination map build
- recombination map population
- gene annotation build and version

If one of these does not match the others, the result may still be technically
generated but biologically misleading.

## Compatibility Matrix

| Layer | What must match | Typical EasyColoc input | Recommended check |
| --- | --- | --- | --- |
| GWAS summary statistics | genome build | `gwas.yaml` dataset `build` | Confirm `hg19` vs `hg38` before harmonization |
| QTL summary statistics | genome build | `qtl.yaml` plus QTL files | Ensure regional queries and phenotype coordinates use the same build |
| LD reference panel | genome build and population | `plink_hg19` or `plink_hg38`, optional `plink_keep` | Match build exactly; prefer ancestry-matched keep/sample files |
| recombination map | genome build and population | `config/global.yaml -> recom` | Use for plotting only, but keep build/pop aligned with the locus panel |
| gene annotation | genome build and version | `config/global.yaml -> gene_anno` | Keep annotation build aligned with the locus coordinates; record version such as GENCODE v40 |
| allele frequency reference | genome build and population | `1kg_af` tables | Match the harmonization build and the intended ancestry when possible |

## Practical Rules

### Build matching

- Do not mix `hg19` GWAS coordinates with `hg38` QTL windows unless the GWAS
  side has been explicitly lifted over and validated.
- The LD panel build must match the coordinates used for clumping and SuSiE LD
  extraction.
- Gene tracks and recombination overlays should be interpreted as annotation
  layers on top of the plotted coordinate system. If the plotted coordinates are
  `hg38`, both layers should also be `hg38`.

### Population matching

- LD panels are population-sensitive and can materially affect fine-mapping and
  credible-set interpretation.
- Recombination maps are less critical than LD for inference, but mismatched
  populations can still make locus context look cleaner or noisier than it
  should.
- For plotting fixtures and demos, population labels should still be explicit so
  the repository does not imply that population choice is arbitrary.

### Version matching

- Gene annotations should record version, not just source.
- A figure that uses `GENCODE v40 hg38` should not be described generically as
  "GENCODE" in documentation.
- If a demo or fixture uses a reduced synthetic annotation derived for testing,
  that should be stated clearly instead of implying it is a real production
  annotation release.

## EasyColoc Smoke Fixture

The plotting smoke fixture is intentionally small but scientifically aligned:

| Component | Fixture value | Why it was chosen |
| --- | --- | --- |
| coordinate system | `hg38` | matches the plotting window encoded in the synthetic phenotype label |
| chromosome | `chr1` | matches the synthetic locus and keeps the fixture small |
| phenotype anchor | `ENST00000654683.1|CCDC30|chr1:42482663-42484158|+` | ensures the transcript window and gene track agree |
| recombination population | `CHB` | explicitly encoded in the recombination fixture path and filename |
| gene annotation version | synthetic test annotation for `hg38` | avoids implying a full production annotation while preserving coordinate correctness |

Relevant fixture files:

- `tests/fixtures/annotation/smoke_hg38_chr1.gtf`
- `tests/fixtures/recomb/hg38/CHB/CHB_recombination_map_hapmap_format_hg38_chr_1.txt`

## Recommended Reporting Practice

For serious analyses, record the following in either the run manifest, the
output report, or both:

- GWAS build
- QTL build
- LD panel build
- LD panel population
- recombination map population
- gene annotation source and version

If these are not recorded, later readers cannot reliably judge whether an
observed difference is biological or reference-driven.
