# Reference Data And Speed

For build, population, and annotation-version matching rules, see
[REFERENCE_COMPATIBILITY.md](REFERENCE_COMPATIBILITY.md).

## Why This Matters

For EasyColoc, the practical adoption bottlenecks are:

1. Users do not know the minimum reference data they actually need.
2. Users do not want to spend hours running unnecessary QTL scans.

This repo now addresses both with:

- `./easycoloc refs`: machine-readable human-facing reference inventory
- `./easycoloc bootstrap-refs`: materialize shared references into a standard project layout
- `analysis.prefilter_sig_pairs: true`: faster QTL phenotype pruning before full coloc

## Quick Inventory

Run from the repository root:

```bash
./easycoloc refs
./easycoloc refs --include-qtl-files
./easycoloc bootstrap-refs config/reference_sources.template.yaml --rewrite-config
./easycoloc bootstrap-refs --demo ./demo_quickstart --run
./easycoloc bootstrap-refs --setup-1kg ./refs/1kg_phase3_hg19 --build hg19 --pop EAS --chromosomes 22
./easycoloc bootstrap-refs --setup-1kg ./refs/1kg_phase3_hg38 --build hg38 --pop EAS --chromosomes 22
./easycoloc bootstrap-refs --fetch-gtex-meta ./refs/gtex_meta --gtex-eqtl-dir /path/to/gtex/eqtl --gtex-sqtl-dir /path/to/gtex/sqtl
```

This prints, for each reference:

- whether it is required or optional
- which pipeline stage uses it
- expected format
- configured path
- whether it is currently present
- current on-disk size

The bootstrap command uses a source-mapping YAML to symlink or copy shared
reference assets into a project-local `refs/` layout and can optionally rewrite
`config/global.yaml` to point at those local paths.

Additional bootstrap modes:

- `--demo DEST_DIR`
  Creates a self-contained toy project with a chr22 PLINK panel, mock GWAS/QTL
  inputs, tabix-indexed QTL files, and project configs. `--run` executes the
  full EasyColoc pipeline and produces an HTML report.

- `--setup-1kg DEST_DIR --build hg19|hg38 --pop EAS`
  Downloads build-specific 1000 Genomes Phase 3 VCFs from the official FTP with
  resumable transfers, converts them to PLINK using `plink`, and generates a
  population-specific `.sample` keep file.

- `--fetch-gtex-meta DEST_DIR`
  Downloads the GTEx sample attributes file. If `--gtex-eqtl-dir` and/or
  `--gtex-sqtl-dir` are provided, it also scans local GTEx QTL files and builds
  summary CSVs that match EasyColoc's `qtl_info.file` expectations.

## Minimum Core References

### Always required

- `plink_hg38`
  Format: PLINK prefix with `.bed/.bim/.fam`
  Used for: locus clumping and SuSiE LD matrix extraction

- `qtl_summary`
  Format: CSV metadata table
  Used for: locating allPairs/sigPairs tabix files

- `1kg_af_*`
  Format: bgzip-compatible allele frequency table per build and population
  Used for: native GWAS harmonization

### Required when hg19 GWAS are present

- `ref_genome_hg19`
  Format: FASTA/FASTA.GZ

- `dbsnp_hg19`
  Format: dbSNP reference file used by harmonization

### Required when hg38 GWAS are present

- `ref_genome_hg38`
  Format: FASTA/FASTA.GZ

- `dbsnp_hg38`
  Format: dbSNP reference file used by harmonization

## Optional But High Value

- `hash_table_dir`
  Value: rescues rsID-to-position matching when direct overlap is poor

- `plink_keep`
  Value: population-specific LD panels for better SuSiE behavior

- `gene_anno`
  Value: gene tracks in locus plots

- `recombination_map`
  Value: recombination ribbon in plots

- `sigPairs` tabix files
  Value: fastest way to avoid unnecessary allPairs coloc runs

## Performance Defaults

EasyColoc now enables this by default in `config/global.yaml`:

```yaml
analysis:
  prefilter_sig_pairs: true
```

What this does:

- Query `sigPairs` first for the locus window
- If no significant phenotype is present, skip that QTL dataset early
- If significant phenotypes are present, only query/filter relevant phenotypes from `allPairs`

This usually reduces wasted QTL work substantially when:

- QTL panels are large
- phenotype counts are high
- only a small subset of phenotypes are signal-bearing at a locus

## Additional Implemented Speedups

- PLINK `.bim` is now cached in-memory per reference prefix instead of being re-read for each GWAS dataset during locus identification.
- QTL metadata paths in the summary CSV are resolved relative to the summary file location, which makes external project layouts more portable and reduces path-debugging time.

## Recommended User Workflow

1. Run `./easycoloc refs --include-qtl-files`
2. Fill only the references that are truly required for your current GWAS build/population mix
3. Keep `prefilter_sig_pairs: true`
4. Run `./easycoloc doctor`
5. Launch with `./easycoloc run --managed`
