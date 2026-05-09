# EasyColoc Docker Quick Start Guide

## Prerequisites

- Docker >= 20.10
- Docker Compose (optional, for easier management)

## Quick Start

### 1. Build the Image

```bash
cd /path/to/EasyColoc
docker build -t easycoloc:latest .
```

**Build time:** ~15-20 minutes (downloads many dependencies)

**Image size:** ~5-6 GB

### 2. Prepare Your Data

Create a data directory with your input files:

```bash
mkdir -p /path/to/my_analysis/{data,reference,output}
```

### 3. Run Analysis

#### Single GWAS Analysis

```bash
docker run -it \
  -v /path/to/my_analysis/data:/data \
  -v /path/to/my_analysis/output:/results \
  -v /path/to/reference:/ref:ro \
  easycoloc:latest
```

#### Batch Processing

```bash
docker run -it \
  -v /path/to/my_analysis:/work \
  easycoloc:latest \
  bash /opt/EasyColoc/tools/batch_run.sh /work/config.yml
```

### 4. Interactive Mode (Debug/Development)

```bash
docker run -it \
  -v /path/to/my_analysis/data:/data \
  easycoloc:latest bash
```

## Volume Mounts Explained

| Host Path | Container Path | Purpose |
|-----------|----------------|---------|
| `/path/to/data` | `/data` | Input GWAS/QTL files |
| `/path/to/results` | `/results` | Analysis output |
| `/path/to/reference` | `/ref` | Reference genomes (read-only) |

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `N_CORES` | `4` | Number of CPU cores for parallel processing |
| `MEMORY_LIMIT` | `8g` | Memory limit for R processes |
| `OUTPUT_DIR` | `/results` | Output directory |

Example with custom settings:

```bash
docker run -it \
  -e N_CORES=8 \
  -e MEMORY_LIMIT=16g \
  -v /data:/data \
  easycoloc:latest
```

## Pre-built Reference Data

For production use, prepare reference data on the host:

```bash
# Directory structure
/reference/
├── 1kg_hg38/           # 1000 Genomes PLINK reference
├── dbsnp/              # dbSNP files
├── gencode.v40.gtf     # Gene annotations
├── recomb/             # Recombination maps
└── liftover/           # UCSC chain files
```

Then mount:

```bash
docker run -v /reference:/ref:ro easycoloc:latest
```

## Docker Compose (Optional)

For complex setups, use `docker-compose.yml`:

```yaml
version: '3.8'

services:
  easycoloc:
    image: easycoloc:latest
    volumes:
      - ./data:/data
      - ./results:/results
      - ./ref:/ref:ro
    environment:
      - N_CORES=8
    deploy:
      resources:
        limits:
          memory: 16G
```

Run with:

```bash
docker-compose up
```

## Troubleshooting

### Permission Issues

If you get permission errors on mounted volumes:

```bash
# Run as same UID/GID as host user
docker run -u $(id -u):$(id -g) \
  -v /data:/data \
  easycoloc:latest
```

### Memory Issues

Limit container memory:

```bash
docker run -m 16g easycoloc:latest
```

### Network Issues

If downloads fail during build, use a mirror:

```dockerfile
# Add to Dockerfile before apt-get/pip install
RUN echo "Acquire::http::Pipeline-Depth 0;" > /etc/apt/apt.conf.d/99custom
```

## Exporting Results

After analysis completes:

```bash
# Copy results from container
docker cp container_id:/results ./my_results

# Or mount volume as shown above
```

## Customization

### Add Custom Reference Data

```dockerfile
FROM easycoloc:latest

COPY my_reference.tar.gz /opt/
RUN tar -xzf /opt/my_reference.tar.gz -C /opt/
```

### Install Additional Packages

```dockerfile
FROM easycoloc:latest

# R packages
RUN R -e "install.packages('your_package', repos='https://cloud.r-project.org')"

# Python packages
RUN pip install your_python_package
```

## Version History

| Tag | R Version | Description |
|-----|-----------|-------------|
| `latest` | 4.3.2 | Latest stable release |
| `1.0` | 4.3.2 | Initial release |
| `dev` | 4.4.0 | Development version |

## Support

- GitHub Issues: https://github.com/cupcake777/EasyColoc/issues
- Documentation: https://github.com/cupcake777/EasyColoc/wiki
