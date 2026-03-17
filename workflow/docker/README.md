# Docker Containers

This directory contains Dockerfiles for tools that are difficult to install via conda or require specific legacy environments.

## Available Containers

### `mafTools/` - MAF Alignment Tools
**Purpose**: Python 2.7-dependent mafTools utilities  
**Base Image**: `python:2.7-slim`  
**Used in**: Step 1 (Extract Ancestor)

**Installed Tools**:
- `mafTools` - MAF alignment manipulation (Python 2.7 only)
- `biopython` (1.76) - Compatible with Python 2.7
- Legacy alignment utilities

**Usage**:
```bash
docker run -v $(pwd):/data juliahoglund/maftools:latest \
  maf_sort.py input.maf > output.maf
```

mafTools was written for Python 2.7 and is not compatible with Python 3.x. Rather than maintaining dual Python installations, it is recommended to use fully containerized versions.

---

### `GERP/` - GERP Conservation Scoring
**Purpose**: GERP (Genomic Evolutionary Rate Profiling) tools  
**Base Image**: `biocontainers/gerp:latest`  
**Used in**: Step 5 (Annotate Variants - Conservation)

**Installed Tools**:
- `gerpcol` - Compute GERP scores from alignments
- `gerpelem` - Element-based conservation

**Usage**:
```bash
docker run -v $(pwd):/data juliahoglund/gerp:latest \
  gerpcol -f alignment.maf -e rates.txt
```

GERP compilation can be problematic on modern systems due to legacy C++ code and library dependencies.

---

### `sift4g/` - SIFT4G Deleteriousness Scoring (to be implemented!)
**Purpose**: SIFT 4G algorithm for variant deleteriousness prediction  
**Base Image**: `biocontainers/sift4g:latest`  
**Used in**: Step 5 (Annotate Variants)

**Installed Tools**:
- `sift4g` - Fast SIFT implementation
- Database utilities

**Usage**:
```bash
docker run -v $(pwd):/data juliahoglund/sift4g:latest \
  sift4g -q query.fasta -d database.fa
```
 
SIFT4G requires specific database formats and compilation flags. Containerization ensures reproducibility.

---

### `renv/` - Complete R Environment
**Purpose**: Full R statistical environment with all packages  
**Base Image**: `rocker/r-ver:4.3`  
**Used in**: Step 4 (Summary Report)

**Installed R Packages**:
- `rmarkdown` - Report generation
- `ggplot2`, `ggtree` - Visualization
- `tidyverse` - Data manipulation
- `BiocManager` - Bioconductor packages
- `GenomicRanges`, `Biostrings` - Genomic utilities
- `DT`, `plotly` - Interactive visualizations

**Usage**:
```bash
docker run -v $(pwd):/data juliahoglund/report-renv:latest \
  Rscript -e "rmarkdown::render('/data/report.Rmd')"
```

R package installation can be fragile, especially for Bioconductor packages with compiled code. The container provides a pre-built, tested environment.

---

## Using Containers in Snakemake

### Docker
```python
# In Snakemake rule
container: "docker://juliahoglund/maftools:latest"
```

### Singularity (HPC Clusters)
```bash
# Snakemake automatically converts Docker → Singularity
snakemake --use-singularity
```

## Pre-built Container Registry

All containers are available on Docker Hub:

| Container | Docker Hub | Size |
|-----------|------------|------|
| mafTools | `juliahoglund/maftools:latest` | ~200 MB |
| GERP | `juliahoglund/gerp:latest` | ~150 MB |
| SIFT4G | `juliahoglund/sift4g:latest` | ~180 MB |
| R environment | `juliahoglund/report-renv:latest` | ~1.2 GB |

## System Requirements

### For Building
- **Docker Desktop** (macOS/Windows) or **Docker Engine** (Linux)
- **Disk space**: 5+ GB for all containers
- **RAM**: 4+ GB recommended

### For Running
- **Docker**: Docker Desktop or Engine
- **Singularity/Apptainer**: v3.0+ (for HPC clusters)

## Future Improvements

- [ ] Add container versioning strategy
- [ ] Implement automated testing for containers
- [ ] Add GPU support for model training containers
