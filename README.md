# EasyColoc

Fast and easy-to-use colocalization analysis pipeline for molecular QTL studies.

## Overview

EasyColoc is an R-based pipeline that performs colocalization analysis between GWAS summary statistics and molecular QTL (eQTL, pQTL, sQTL, etc.) datasets. It features automatic SNP matching, build conversion, multiple colocalization methods, and publication-ready visualization.

## Features

- **Multiple Matching Strategies**: rsID direct match, chromosome-position conversion, position+allele matching
- **Automatic Build Conversion**: Built-in harmonization with gwaslab (hg19 ↔ hg38)
- **Multiple Coloc Methods**: Traditional ABF method and SuSiE-based fine-mapping
- **Interactive Visualization**: LocusZoom-style plots with recombination tracks and gene annotations
- **Parallel Processing**: Efficient multi-core processing for large datasets

## Data Flow Architecture

EasyColoc uses a **two-pipeline architecture** to handle differences between GWAS and QTL data formats:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         GWAS PIPELINE (Left Side)                           │
└─────────────────────────────────────────────────────────────────────────────┘
   Raw GWAS Summary Statistics
   (May have: missing rsID, hg19 coordinates, strand ambiguity)
           │
           ▼
   ┌───────────────────────────┐
   │  gwaslab Harmonization   │  ← rsID annotation, liftOver hg19→hg38,
   │  (Python preprocessing)   │    strand flip, allele standardization
   └───────────────────────────┘
           │
           ▼
   Clean GWAS Data (hg38)
   ✓ rsID annotated (99.98% coverage)
   ✓ hg38 coordinates
   ✓ Standardized alleles (A/C/G/T)
   Format: rsID-based (e.g., "rs123456")
           │
           │
           ▼
   ┌─────────────────────────────────────────────────────────────────────────┐
   │                    THREE-TIER MATCHING STRATEGY                          │
   ├─────────────────────────────────────────────────────────────────────────┤
   │                                                                           │
   │  [Strategy 1] rsID Direct Match                                          │
   │  ├─ Try: Match GWAS rsID ↔ QTL rsID                                      │
   │  └─ Success: ~5% of cases (when QTL provides rsID)                       │
   │                                                                           │
   │  [Strategy 2] Hash Table Conversion  ← WHY HASH_TABLE IS NEEDED          │
   │  ├─ Problem: QTL uses chr:pos (e.g., "1_12345"), GWAS uses rsID          │
   │  ├─ Solution: rsID → chr:pos conversion via dbSNP hash tables             │
   │  ├─ Process: Load chromosome-specific RDS → lookup rsID → get position   │
   │  └─ Success: ~90% of cases (standard QTL sources like GTEx)              │
   │                                                                           │
   │  [Strategy 3] Position + Allele Match                                    │
   │  ├─ Fallback: Match by CHR_POS_EA_NEA variant identifiers                │
   │  ├─ Handles: Novel variants missing from dbSNP                            │
   │  └─ Success: ~5% of cases (edge cases + missing rsIDs)                   │
   │                                                                           │
   └─────────────────────────────────────────────────────────────────────────┘
           │
           ▼
   Merged Dataset (GWAS + QTL)
   → Colocalization Analysis (coloc ABF / SuSiE)
   → Visualization (LocusZoom plots)


┌─────────────────────────────────────────────────────────────────────────────┐
│                          QTL PIPELINE (Right Side)                          │
└─────────────────────────────────────────────────────────────────────────────┘
   Raw QTL Data (Tabix-indexed)
   Source: GTEx, eQTL Catalogue, BLUEPRINT, etc.
   Format: chr:pos-based (e.g., "1_12345" or "1_12345_A_G")
           │
           ▼
   ┌───────────────────────────┐
   │   Direct Loading          │  ← NO gwaslab processing
   │   (Tabix query)           │    (already in standard format)
   └───────────────────────────┘
           │
           ▼
   Clean QTL Data (hg38)
   ✓ Chromosome:Position format
   ✓ hg38 coordinates (by source)
   ✓ Pre-computed effect sizes
   Format: chr:pos-based (NO rsID)
           │
           └─────────────────────────────────────┐
                                                 │
                                                 ▼
                                    (Merges with GWAS via
                                     Three-Tier Strategy)
```

### Key Design Decisions

1. **Why QTL bypasses gwaslab**: QTL consortia (GTEx, eQTL Catalogue) provide pre-harmonized data in standardized chr:pos format. Re-processing would:
   - Waste computational resources (redundant liftOver)
   - Risk introducing errors (unnecessary transformations)
   - Lose consortium-specific quality control

2. **Why hash_table is NOT redundant**: 
   - **Problem**: GWAS uses rsID ("rs123456"), QTL uses chr:pos ("1_12345")
   - **Solution**: hash_table bridges the format gap via dbSNP lookups
   - **When triggered**: Automatically when Strategy 1 (rsID direct match) fails to find ≥30 SNPs
   - **Performance**: Lazy loading per chromosome (~0.5s lookup time)

3. **Graceful degradation**: If hash_table is unavailable or incomplete:
   - Falls back to Strategy 3 (position + allele matching)
   - Ensures pipeline never fails due to missing reference data

## Installation

### Environment Setup (Recommended)

EasyColoc provides a complete conda environment configuration for easy installation.

```bash
# Install using micromamba (recommended)
micromamba create -f tools/easycoloc.yml
micromamba activate easycoloc

# Or using conda
conda env create -f tools/easycoloc.yml
conda activate easycoloc
```

### Manual Installation

#### System Requirements

- R >= 4.2.0
- Python >= 3.9
- PLINK >= 1.9
- Tabix
- BCFtools

#### R Dependencies

```r
install.packages(c(
  "data.table", "dplyr", "yaml", "parallel", "tools",
  "glue", "coloc", "susieR", "ggplot2", "ggpubr",
  "GenomicRanges", "vroom", "forcats", "purrr",
  "ggrepel", "rtracklayer", "clusterProfiler",
  "org.Hs.eg.db", "jsonlite", "tidyr",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "GenomicFeatures", "AnnotationDbi"
))
```

#### Python Dependencies

```bash
pip install gwaslab==3.6.3 adjusttext biopython gtfparse liftover \
  matplotlib numpy pandas pyensembl pysam scikit-allel \
  scipy seaborn h5py pyarrow polars
```

### Quick Setup

```bash
git clone https://github.com/cupcake777/EasyColoc.git
cd EasyColoc
micromamba create -f tools/easycoloc.yml
micromamba activate easycoloc
```

## Data Preparation

Before running EasyColoc, you need to prepare reference data files. All paths should be configured in `config/global.yaml`.

### Required Reference Files

| Resource | Description | Download |
|----------|-------------|----------|
| **1000 Genomes** | PLINK format genotype reference (hg38) | https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg |
| **dbSNP** | SNP position mapping (hg38) | https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/ |
| **Recombination Maps** | Genetic recombination rates | https://zenodo.org/records/11437540 |
| **Gene Annotation** | GENCODE GTF file | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/ |

### Build Conversion Reference

If your GWAS data is in hg19, you'll need liftOver chain files:
- hg19 → hg38: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
- hg38 → hg19: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

### Setting Up SNP Hash Tables

EasyColoc uses pre-computed hash tables for fast rsID ↔ position conversion.

1. **Download dbSNP BED files**: Visit [NCBI dbSNP FTP](https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/) and download all `bed_chr_*.bed.gz` files to a local directory.

2. **Convert BED directly to RDS**:
```bash
Rscript tools/bed2rds.r -i /path/to/bed_files -o /path/to/hash_tables
```

**Legacy method** (BED → JSON → RDS):

```bash
python tools/dbsnp_hash_table_maker.py
Rscript tools/hash2rds.r
```

### Configuration Example

```yaml
# config/global.yaml
output_dir: "/path/to/results"
use_parallel: true
n_cores: 8

# Reference data paths
plink_hg38: "/path/to/1kg_hg38/1000g"
plink_hg19: "/path/to/1kg_hg19/1000g"
gene_anno: "/path/to/gencode.v43.gtf"
recomb_map: "/path/to/recomb_rate/CHB/CHB"

# Hash table for rsID conversion
hash_table_dir: "/path/to/snp_hash_tables/"
```

### Test Data

Sample test data is available for download to help you get started:

- **Dropbox**: [EasyColoc Test Data](https://www.dropbox.com/scl/fi/r27ura8jr83j6yoo26phn/test_data.tar.xz?rlkey=mpzjcnejmf2rzxnwyfsvbtd1y&st=4x7k63zj&dl=0)
- Contains: Sample GWAS summary stats, QTL data.

```bash
# Download and extract test data
wget https://www.dropbox.com/scl/fi/r27ura8jr83j6yoo26phn/test_data.tar.xz?rlkey=mpzjcnejmf2rzxnwyfsvbtd1y&st=4x7k63zj&dl=0
tar -xf test_data.tar.xz -C /path/to/EasyColoc/data/
```

## Quick Start

### 1. Configure Data

Edit the configuration files in the `config/` directory:

```yaml
# config/global.yaml
output_dir: "/path/to/results"
use_parallel: true
n_cores: 8
plink_hg38: "/path/to/1kg_hg38"
gene_anno: "/path/to/gencode.gtf"
hash_table_dir: "/path/to/snp_hash_tables/"
```

```yaml
# config/gwas.yaml
datasets:
  - id: "MyGWAS"
    name: "My GWAS Study"
    file: "/path/to/gwas.txt"
    build: "hg19"
    type: "cc"
    prop: 0.5
    columns:
      snp: "SNP"
      chrom: "CHR"
      pos: "POS"
      beta: "BETA"
      se: "SE"
      pval: "P"
      n: "N"
```

```yaml
# config/qtl.yaml
qtl_info:
  file: "/path/to/QTL_summary.csv"
  base_dir: "/path/to/qtl"
  columns:
    id: "Type"
    sample_size: "N"
    all_filename: "allPairsTabixFilename"
```

### 2. Run Pipeline

```bash
Rscript run_coloc.r
```

## Configuration Options

### Global Settings (config/global.yaml)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `output_dir` | Output directory | Required |
| `use_parallel` | Enable parallel processing | true |
| `n_cores` | Number of cores for parallelization | 8 |
| `plink_hg38` | hg38 PLINK reference prefix | Required |
| `gene_anno` | Path to gene annotation (GTF) | Required |
| `hash_table_dir` | rsID conversion tables directory | Optional |

#### Advanced Settings

You can fine-tune analysis and plotting parameters in `config/global.yaml`:

```yaml
# Clumping parameters (for Locus Identification)
clump:
  p1: 5e-8        # Primary P-value threshold (index SNPs)
  p2: 5e-8        # Secondary P-value threshold (clumped SNPs)
  r2: 0.1         # LD r^2 threshold
  kb: 1000        # Clumping window in kb

coloc_settings:
  p1: 1e-4        # Prior probability a SNP is associated with trait 1
  p2: 1e-4        # Prior probability a SNP is associated with trait 2
  p12: 1e-5       # Prior probability a SNP is associated with both traits
  susie_args:     # Arguments for SuSiE fine-mapping
    L: 10         # Maximum number of non-zero effects

harmonization_settings:
  vcf_ref: "/path/to/ref.vcf.gz" # Optional VCF for harmonization
  target_build: "hg38"           # Target genome build

plot_settings:
  style: "locuszoom"   # Plot style: "locuszoom" or "classic"
  width: 10            # Plot width in inches
  height: 6            # Plot height in inches
  res: 300             # DPI resolution
  window: 1000000      # Window size in bp (default 1MB)
  color_scheme: "navy_red" # LD color scheme
```

### GWAS Settings (config/gwas.yaml)

| Parameter | Description | Required |
|-----------|-------------|----------|
| `id` | Dataset identifier | Yes |
| `file` | Path to summary stats | Yes |
| `build` | Genome build (hg19/hg38) | Yes |
| `type` | Study type (cc/quant) | Yes |
| `prop` | Case proportion (cc only) | If type=cc |

## Output

```
output/
├── abf/                    # Coloc results
│   └── *_results.csv
├── susie/                  # SuSiE fine-mapping
│   └── *_susie.csv
└── plots/                  # Visualization plots
    └── *_coloc.pdf
```

### Results Format

| Column | Description |
|--------|-------------|
| Locus | Lead SNP identifier |
| Gene | Gene symbol |
| PP4 | Posterior probability of colocalization (H4) |
| n_snps | Number of SNPs in analysis |

### Interpretation Guide

| PP4 Range | Interpretation |
|-----------|----------------|
| 0 - 0.5 | No evidence of colocalization |
| 0.5 - 0.8 | Suggestive evidence |
| 0.8 - 0.95 | Strong evidence |
| > 0.95 | Very strong evidence |

## Directory Structure

```
EasyColoc/
├── run_coloc.r          # Main pipeline script
├── README.md            # This file
├── LICENSE              # MIT License
├── config/              # Configuration files
│   ├── global.yaml
│   ├── gwas.yaml
│   └── qtl.yaml
├── src/                 # Utility modules
│   ├── utils_analysis.R # Coloc and SuSiE methods
│   ├── utils_format.R   # Data formatting and merging
│   ├── utils_hash.R     # SNP ID conversion
│   ├── utils_helpers.R  # Helper functions
│   └── utils_plot.R     # Visualization
└── tools/               # Utility scripts
    ├── bed2rds.r             # Convert BED directly to RDS
    └── easycoloc.yml        # Conda environment configuration
```

## Citation

EasyColoc is inspired by [ColocQuiaL](https://github.com/bvoightlab/ColocQuiaL) and utilizes [gwaslab](https://github.com/Jia-Xuan-Low/gwaslab) for data harmonization.

If you use EasyColoc in your research, please cite:

```bibtex
@repo{EasyColoc,
  author = {cupcake777},
  title = {EasyColoc: Fast and easy-to-use colocalization analysis},
  url = {https://github.com/cupcake777/EasyColoc},
  year = {2025}
}
```

## License

MIT License - see [LICENSE](LICENSE) file.
