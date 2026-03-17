# Conda Environments

This directory contains conda environment definitions for the omniCADD pipeline. Each environment is specialized for specific workflow steps to minimize dependency conflicts and installation time.

## Environment Files

### `ancestor.yml` - Ancestral Sequence Processing
**Used in**: Step 1 (Extract Ancestor)

**Key Dependencies**:
- `biopython` - Sequence manipulation and parsing
- `pyfaidx` - FASTA indexing and retrieval
- `lz4` - Compression for intermediate files
- `msaconverter` - Multiple sequence alignment format conversion
- `pandas`, `numpy` - Data processing

**Purpose**: Process multi-species alignments to extract ancestral sequences

---

### `annotation.yml` - Variant Annotation Tools
**Used in**: Step 5 (Annotate Variants)

**Key Dependencies**:
- `ensembl-vep` (v112+) - Variant Effect Predictor
- `phast` - phyloFit, phastCons, phyloP conservation tools
- `gerp` - GERP conservation scoring
- `samtools`, `bcftools` - Variant file processing
- `bedtools` - Genomic interval operations
- `sift4g` - SIFT deleteriousness prediction
- `tabix`, `htslib` - File indexing

**Purpose**: Functional and conservation annotation of variants

**Notes**:
- VEP cache may need manual download if auto-install fails
- GERP binaries may require Docker fallback on some systems

---

### `model.yml` - Machine Learning
**Used in**: Step 6 (Combine Annotations), Step 7 (Train Model), Step 8 (Score Variants)

**Key Dependencies**:
- `scikit-learn` - Logistic regression classifier
- `pandas` - DataFrame operations
- `numpy` - Numerical computing
- `scipy` - Sparse matrix operations
- `joblib` - Model serialization
- `matplotlib`, `seaborn` - Visualization

**Purpose**: Feature matrix manipulation, model training, and scoring

**Notes**:
- Supports sparse matrix operations for memory efficiency
- Includes utilities for cross-validation and hyperparameter tuning

---

### `score.yml` - Variant Generation & Filtering
**Used in**: Step 8 (Score Variants)

**Key Dependencies**:
- `bcftools` - VCF manipulation
- `bedtools` - Genomic operations
- `gatk` (3.6) - Variant filtering
- `samtools` - Reference genome processing
- `scikit-learn` - Scoring with trained model
- `tabix` - VCF indexing

**Purpose**: Generate all possible variants and apply CADD scores

**Notes**:
- GATK version locked to 3.6 for compatibility
- Includes VCF utilities for efficient variant generation

---

### `simulation.yml` - Variant Simulation
**Used in**: Step 3 (Simulate Variants)

**Key Dependencies**:
- `biopython` - Sequence manipulation
- `vcftools` - VCF file processing
- `bcftools` - VCF merging and filtering
- `pandas` - Data processing

**Purpose**: Simulate variants enriched for deleterious mutations (positive class)

**Notes**:
- Implements substitution rate models
- Supports both SNP and INDEL (unstable!) simulation

---

### `common.yml` - Utilities
**Used in**: Multiple steps

**Key Dependencies**:
- `bgzip` - Block gzip compression
- `tabix` - Tabix indexing
- `bcftools` - VCF operations
- `vcftools` - VCF utilities
- `samtools` - SAM/BAM/CRAM tools

**Purpose**: Common genomic utilities used across pipeline

**Notes**:
- Lightweight environment for quick operations
- Used for file format conversions and indexing

---

### `r_stats.yml` - R Environment for Reporting
**Used in**: Step 4 (Summary Report)

**Key Dependencies**:
- `r-base` - R statistical environment
- `r-rmarkdown` - Report generation
- `r-ggplot2` - Data visualization
- `r-ggtree` - Phylogenetic tree plotting
- `r-tidyverse` - Data manipulation
- `r-dt` - Interactive tables
- `bioconductor-*` - Genomic R packages

**Purpose**: Generate interactive HTML summary reports

**Notes**:
- Can be replaced with Docker container if installation issues occur
- Requires X11 libraries for some plotting functions

---

## Requirements

### System Requirements
- **Conda** or **Mamba** (mamba recommended for faster installation)
- **Disk space**: ~20 GB for all environments combined
- **RAM**: 4+ GB for environment solving

### Channel Configuration
Ensure proper channel priority in `~/.condarc`:

```yaml
channels:
  - conda-forge
  - bioconda
  - defaults
channel_priority: strict
```

## Environment Isolation

Each environment is isolated to prevent dependency conflicts:

- **VEP** and **PHAST** separated to avoid library conflicts
- **Python 2.7** tools (mafTools) run in Docker
- **R environment** isolated from Python tools

## Future Improvements

- [ ] Add GPU-accelerated variants for model training
- [ ] Add pre-built containers for all environments
- [ ] Update to latest VEP and tool versions
