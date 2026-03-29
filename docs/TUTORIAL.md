# Tutorial

This guide walks from a toy demo to a real analysis setup.

## 1. Environment Setup

### Required tools

EasyColoc expects these tools to be available locally:

- `Rscript`
- `plink`
- `bgzip`
- `tabix`

For the full GWAS harmonization path, you also want:

- a working `gwaslab` environment
- reference FASTA/dbSNP resources

Recommended first checks:

```bash
./easycoloc refs
./easycoloc doctor
```

### Optional quick validation

```bash
./easycoloc smoke
```

## 2. Demo In Under Two Minutes

Create a fully self-contained project and run it immediately:

```bash
./easycoloc bootstrap-refs --demo ./demo_quickstart --run
```

After it finishes, inspect:

- `demo_quickstart/results/coloc_report.html`
- `demo_quickstart/results/all_colocalization_results.csv`
- `demo_quickstart/results/all_susie_results.csv`

This is the fastest way to confirm that your local EasyColoc installation,
`plink`, `bgzip`, and `tabix` are working together.

## 3. Build A Real 1000 Genomes Reference

### hg19 example

```bash
./easycoloc bootstrap-refs \
  --setup-1kg ./refs/1kg_phase3_hg19 \
  --build hg19 \
  --pop EAS \
  --chromosomes 1-22
```

### hg38 example

```bash
./easycoloc bootstrap-refs \
  --setup-1kg ./refs/1kg_phase3_hg38 \
  --build hg38 \
  --pop EAS \
  --chromosomes 1-22
```

What this does:

- downloads build-specific 1000 Genomes VCFs with resumable transfer
- converts them to PLINK
- creates a population-specific `.sample` keep file

If you want EasyColoc to wire the resulting panel into your config automatically,
run the same command with `--rewrite-config` and pass your config paths.

## 4. Prepare GTEx Metadata

If you already have local GTEx QTL files:

```bash
./easycoloc bootstrap-refs \
  --fetch-gtex-meta ./refs/gtex_meta \
  --gtex-eqtl-dir /path/to/gtex/eqtl \
  --gtex-sqtl-dir /path/to/gtex/sqtl
```

This writes:

- `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`
- `GTEx_v8_eQTL_summary.csv`
- `GTEx_v8_sQTL_summary.csv`
- `qtl_gtex_generated.yaml`

Use the generated YAML as a starting point for GTEx-based analyses.

## 5. Prepare Your Configs

The minimum files to check are:

- `config/global.yaml`
- `config/gwas.yaml`
- `config/qtl.yaml`

### `config/global.yaml`

This is where you define:

- output and temp directories
- build-aware PLINK references: `plink_hg19`, `plink_hg38`
- `plink_keep`
- FASTA / dbSNP / 1KG AF resources
- coloc thresholds
- plotting and runtime behavior

### `config/gwas.yaml`

This defines one or more GWAS datasets:

- input path
- build
- population
- trait type
- sample size / case fraction
- column mapping

### `config/qtl.yaml`

This defines:

- QTL build
- QTL summary CSV
- which columns point to `allPairs` and `sigPairs`
- QTL field mapping for downstream parsing

Common fields to confirm:

- `plink_hg19` / `plink_hg38`
- `plink_keep`
- `1kg_af`
- `dbsnp_hg19` / `dbsnp_hg38`
- `qtl_info.file`

Run:

```bash
./easycoloc refs --include-qtl-files
./easycoloc doctor
```

## 6. Full Pipeline Run

For normal work, prefer the managed wrapper:

```bash
./easycoloc run --managed
```

Or run against explicit config files:

```bash
./easycoloc run --managed \
  --global /path/to/global.yaml \
  --gwas /path/to/gwas.yaml \
  --qtl /path/to/qtl.yaml
```

What happens during the run:

1. EasyColoc harmonizes GWAS input when required.
2. Significant loci are discovered with PLINK clumping.
3. QTL datasets are prefiltered with `sigPairs` when available.
4. Variant matching is performed by rsID, hash rescue, or allele-aware position matching.
5. ABF coloc runs first; SuSiE follows when signal and LD support it.
6. Plots, RDS bundles, merged summaries, and an HTML report are written.

## 7. Inspect Progress And Results

During or after the run:

```bash
./easycoloc status /path/to/output_dir
./easycoloc monitor /path/to/output_dir
./easycoloc check /path/to/output_dir
./easycoloc manifest /path/to/output_dir
```

Main result files:

- `all_colocalization_results.csv`
- `significant_colocalizations_PP4_*.csv`
- `all_susie_results.csv`
- `coloc_report.html`

## 8. Regenerate Plots Without Re-running Coloc

If RDS bundles already exist:

```bash
Rscript tools/rerun_plots.R
```

This regenerates locus plots from saved RDS objects without rerunning the whole
pipeline.

## 9. Common Failure Modes

### `plink` missing

Bootstrap and 1KG setup require `plink` in `$PATH`. EasyColoc now fails with an
explicit install hint when it is missing.

### QTL files present but no output

Check that:

- `qtl_info.file` points to the correct summary CSV
- `allPairsTabixFilename` and `sigPairsTabixFilename` resolve correctly
- `sigPairs` are not filtering everything unexpectedly at the locus of interest

### Build mismatch

Make sure the QTL build matches the PLINK reference build selected in
`config/global.yaml`:

- `plink_hg19` for `qtl_info.build: hg19`
- `plink_hg38` for `qtl_info.build: hg38`
