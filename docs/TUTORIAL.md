# Tutorial

This guide is organized by common usage scenarios so you can jump directly to
the path you need.

## Scenario 1. Open An Existing Results Directory

If you already have a completed `results/` directory from EasyColoc, start here:

```bash
./easycoloc report-web /path/to/results
```

This command prepares report payload data and starts a local report server for
interactive browsing.

## Scenario 2. First-Time Environment Setup

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
micromamba create -f environment.yml
micromamba activate easycoloc
./easycoloc smoke
```

## Scenario 3. Demo In Under Two Minutes

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

## Scenario 4. Choose Your Starting Mode

| Situation | Best next step |
| --- | --- |
| You only want to validate the repo | Run `./easycoloc smoke` |
| You need a portable starter project | Run `./easycoloc init <dir>` |
| You already have reference assets | Edit the portable defaults in `config/*.yaml` or create private overrides in `config/local/` |
| You need local 1KG or GTEx support files | Use `bootstrap-refs` as shown below |

## Scenario 5. Build A Real 1000 Genomes Reference

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

## Scenario 6. Prepare GTEx Metadata

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

## Scenario 7. Prepare Your Configs

The minimum files to review are:

- `config/global.yaml`
- `config/gwas.yaml`
- `config/qtl.yaml`
- `config/README.md`

| File | Main responsibility |
| --- | --- |
| `config/global.yaml` | output paths, references, coloc thresholds, plotting, runtime behavior |
| `config/gwas.yaml` | GWAS datasets, build, population, trait type, sample size, column mapping |
| `config/qtl.yaml` | QTL build, metadata table, allPairs/sigPairs columns, downstream field mapping |

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

If you do not want to commit machine-specific paths, put private copies in
`config/local/` and pass them with `--global`, `--gwas`, and `--qtl`.

## Scenario 8. Full Pipeline Run

For routine work, prefer the managed wrapper:

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

## Scenario 9. Inspect Progress And Results

During or after the run:

```bash
./easycoloc status /path/to/output_dir
./easycoloc monitor /path/to/output_dir
./easycoloc watch /path/to/output_dir 60 logs/monitor/current_run.log
./easycoloc check /path/to/output_dir
./easycoloc manifest /path/to/output_dir
```

`watch` is useful when the results directory is outside the current project
tree or should remain read-only during auditing. It writes snapshots to a local
log file instead of mutating the target output directory.

| File | Why you care |
| --- | --- |
| `all_colocalization_results.csv` | complete merged coloc table |
| `significant_colocalizations_PP4_*.csv` | thresholded candidate hits |
| `all_susie_results.csv` | merged fine-mapping output |
| `coloc_report.html` | quick visual review for handoff or audit |

## Scenario 10. Regenerate Plots Without Re-running Coloc

If RDS bundles already exist:

```bash
Rscript tools/rerun_plots.R
```

This regenerates locus plots from saved RDS objects without rerunning the whole
pipeline.

## Scenario 11. Common Failure Modes

| Symptom | What to check |
| --- | --- |
| `plink` missing | Ensure `plink` is on `$PATH`; bootstrap and clumping require it |
| QTL files present but no output | Confirm `qtl_info.file`, `allPairsTabixFilename`, and `sigPairsTabixFilename`, and verify `sigPairs` is not filtering every locus |
| Build mismatch | Match `qtl_info.build: hg19` to `plink_hg19`, and `qtl_info.build: hg38` to `plink_hg38` |
