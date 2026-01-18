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

```bash
# Download dbSNP BED files (choose one build)
# hg38:
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/bed_chr_*.bed.gz

# Convert BED directly to RDS (single step)
Rscript tools/bed2rds.r -i /path/to/bed_files -o /path/to/hash_tables

# Or with options
Rscript tools/bed2rds.r --input /path/to/bed_files \
                        --output /path/to/hash_tables \
                        --build hg38
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
plink_hg38: "/path/to/1kg_hg38/1000g.pgen"
plink_hg19: "/path/to/1kg_hg19/1000g.pgen"
gene_anno: "/path/to/gencode.v43.gtf"
recomb_map: "/path/to/recomb_rate/"

# Hash table for rsID conversion
hash_table_dir: "/path/to/snp_hash_tables/"
```

### Test Data

Sample test data is available for download to help you get started:

- **Dropbox**: [EasyColoc Test Data](https://www.dropbox.com/sh/xxxxx)
- Contains: Sample GWAS summary stats, QTL data, and configuration files

```bash
# Download and extract test data
wget https://www.dropbox.com/sh/xxxxx/test_data.tar.gz
tar -xzf test_data.tar.gz -C /path/to/easycoloc/
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
| `sig_threshold` | PP4 threshold for plotting | 0.8 |
| `use_parallel` | Enable parallel processing | true |
| `plink_hg38` | hg38 PLINK reference | Required |
| `hash_table_dir` | rsID conversion tables | Optional |

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
    ├── bed2rds.r             # Convert BED directly to RDS (recommended)
    ├── hash2rds.r           # Convert JSON to RDS (legacy)
    ├── dbsnp_hash_table_maker.py  # Create hash tables from BED files (legacy)
    └── easycoloc.yml        # Conda environment configuration
```

## Citation

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
