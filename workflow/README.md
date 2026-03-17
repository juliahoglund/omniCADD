# omniCADD Workflow

This directory contains all components of the omniCADD pipeline: Snakemake rules, processing scripts, conda environments, and Docker containers.

## Workflow Structure

```
workflow/
├── rules/          # Snakemake rule files (one per workflow step)
├── scripts/        # Python/Perl/R scripts organized by step
├── envs/           # Conda environment YAML files
├── docker/         # Dockerfiles for problematic dependencies
└── README.md       # This file
```

## Pipeline Steps

The omniCADD workflow consists of 8 steps:

### Step 1: Extract Ancestral Sequences
**Rule**: `1_extract_ancestor.smk`  
**Purpose**: Reconstruct most recent common ancestor (MRCA) sequences from multi-species alignments

**Process**:
- Clean MAF alignment files (remove ambiguous nucleotides, gaps)
- Identify ancestral node between species of interest and comparison species
- Deduplicate sequences, reorder species, correct strand orientation
- Sort alignment blocks and extract ancestor sequences

**Key Scripts**: `step_1_extract_ancestor/*.py`  
**Inputs**: EPO alignments (MAF format), phylogenetic tree  
**Outputs**: `results/ancestral_seq/{ancestor_name}/chr*.fa`

---

### Step 2: Derive Variants
**Rule**: `2_derive_variants.smk`  
**Purpose**: Identify derived variants from population VCF data

**Process**:
1. Extract high-frequency variants (>90% AF) from population VCF
2. Apply 5-case logic comparing ancestor, reference, and population alleles
3. Filter for singleton SNPs (remove clustered variants)

**Key Scripts**: `step_2_derive_variants/*.py`  
**Inputs**: Population VCF, reference genome, ancestral sequences  
**Outputs**: `results/derived_variants/singletons/chr*.vcf`  
**Note**: Assumes population-level VCF with allele frequencies

---

### Step 3: Simulate Variants
**Rule**: `3_simulate_variants.smk`  
**Purpose**: Generate simulated variants enriched for deleterious mutations for model training

**Process**:
1. Compute mutation parameters (substitution rates from ancestor ↔ reference)
2. Scale parameters by overestimation factor (default 5x)
3. Simulate SNPs and optionally INDELs per chromosome
4. Filter variants on gaps or missing ancestral sequence

**Key Scripts**: `step_3_simulate_variants/*.py`  
**Inputs**: Reference genome, ancestral sequences  
**Outputs**: `results/simulated_variants/trimmed_{snps|indels}/chr*.vcf`  
**Algorithm**: Based on Kircher et al. 2014 original implementation

---

### Step 4: Summary Report
**Rule**: `4_summary_report.smk`  
**Purpose**: Generate HTML visualization of genome and variant properties

**Outputs**:
- Variant statistics (SNP counts, transition/transversion ratios)
- Ideogram-style genome graph
- Phylogenetic tree visualization
- Optional CDS coverage annotations (if GFF3/GTF provided)

**Key Scripts**: `step_4_summary_report/*.Rmd`  
**Tech Stack**: R with rmarkdown, ggplot2, ggtree  
**Outputs**: `output/{species_comparison}.html`  
**Note**: Can run in Docker container if R packages unavailable locally

---

### Step 5: Annotate Variants
**Rules**: `5_annotate_vars.smk`  
**Purpose**: Apply functional and conservation annotations to variants

**Annotation Types**:

**A. Functional Annotation (VEP/SNPEff)**:
- Variant effect prediction (missense, nonsense, etc.)
- SIFT deleteriousness scores (when available, to be implemented)
- Grantham matrix amino acid substitution severity

**B. Conservation - PHAST**:
- `phyloFit`: Train phylogenetic substitution models
- `phastCons`: Conservation scores
- `phyloP`: P-values for accelerated/conserved evolution

**C. Conservation - GERP**:
- `gerpcol`: Compute expected vs observed substitution rates
- Map GERP rates to genomic coordinates

**Key Scripts**: `step_5_annotate_variants/*.py`  
**Inputs**: Derived/simulated variants, reference, alignments  
**Outputs**: `results/annotations/{vep|phast|gerp}/*.tsv`

---

### Step 6: Combine Annotations
**Rule**: `6_combine_annotations.smk`  
**Purpose**: Merge all annotations into unified feature matrices

**Process**:
1. Split alignments into manageable chunks (~10k positions)
2. Convert MAF → FASTA format
3. Prune alignment columns with excessive gaps in reference genome (to maintain correct bp positioning in downstream steps)
4. Run GERP and PHAST on chunks in parallel
5. Merge VEP + GERP + PHAST outputs
6. Create sparse numpy NPZ datasets with feature matrices
7. Calculate imputed values (configurable) for missing values
8. Analyze feature relevance

**Key Scripts**: `step_6_combine_annotations/*.py`  
**Outputs**:
- `results/dataset/{derived|simulated}/chr*.npz` (sparse matrices)
- `results/dataset/imputation_dict.txt` (feature means)
- `results/figures/column_analysis/relevance.tsv`

---

### Step 7: Train & Test Model
**Rule**: `7_train_test_model.smk`  
**Purpose**: Build logistic regression classifier to predict deleteriousness

**Process**:
1. Split combined dataset into K folds
2. K-fold cross-validation with grid search:
3. Train final model on all data with best hyperparameters
4. Export model, scaler, and feature coefficients

**Key Scripts**: `step_7_train_test_model/*.py`  
**Algorithm**: scikit-learn LogisticRegression with L2 regularization  
**Outputs**:
- `results/models/final_model.pkl`
- `results/models/feature_scaler.pkl`
- `results/models/feature_coefficients.tsv`

---

### Step 8: Score Variants
**Rule**: `8_score_variants.smk`  
**Purpose**: Generate genome-wide PHRED-scaled CADD scores

**Process**:
1. Generate all possible SNPs per chromosome in blocks (~10k positions)
2. Annotate all variants with VEP
3. Intersect with conservation scores (PHAST, GERP)
4. Prepare feature matrices for prediction
5. Score variants using trained model (deleteriousness probability)
6. Convert to PHRED scale: `PHRED = -10 * log10(1 - rank/total)`
7. Merge blocks, sort, compress (bgzip), and index (tabix)

**Key Scripts**: `step_8_score_variants/*.py`  
**Outputs**:
- `results/cadd_scores/chr*.tsv.gz` (tabix-indexed)
- `results/cadd_scores/scoring_summary.txt`

---

## Data Flow

```
External Resources → Step 1: Extract Ancestor
                           ↓
        ┌──────────────────┴──────────────────┐
        ↓                                      ↓
Step 2: Derive Variants              Step 3: Simulate Variants
   (Real/Deleterious)                   (Neutral/Benign)
        ↓                                      ↓
        └──────────────────┬───────────────────┘
                           ↓
                  Step 5: Annotate
                 (VEP, PHAST, GERP)
                           ↓
              Step 6: Combine Annotations
                    (Feature Matrices)
                           ↓
             Step 7: Train Model (ML Classifier)
                           ↓
              Step 8: Score All Variants
                           ↓
                    CADD Scores
```

Step 4 (Summary Report) runs independently after Steps 1-3.

## Conda Environments

See [envs/README.md](envs/README.md) for detailed environment documentation.

| Environment | Key Tools | Purpose |
|-------------|-----------|---------|
| `ancestor.yml` | biopython, pyfaidx, msaconverter | Alignment processing |
| `annotation.yml` | ensembl-vep, phast, gerp | Variant annotation |
| `model.yml` | scikit-learn, pandas, numpy | Machine learning |
| `score.yml` | bcftools, bedtools, gatk | Variant generation |
| `simulation.yml` | biopython, vcftools | Variant simulation |
| `common.yml` | bgzip, tabix, bcftools | Utilities |
| `r_stats.yml` | R, rmarkdown, ggplot2 | Reporting |

## Docker Containers

Some dependencies require legacy software or complex builds. See [docker/README.md](docker/README.md).

| Container | Purpose | Base Image |
|-----------|---------|------------|
| `mafTools` | Python 2.7 mafTools | Python 2.7 |
| `GERP` | GERP conservation | biocontainers |
| `sift4g` | SIFT4G scoring | biocontainers |
| `renv` | Complete R environment | rocker/r-ver |

## Current Limitations

- ⚠️ **Alignment-free mode**: Rules not yet revised for alignment-free branch
- ⚠️ **VEP cache**: Auto-installation may fail; manual download may be required
- ⚠️ **Memory**: Large chromosomes may require 64+ GB RAM for annotation
- ⚠️ **Python 2.7**: mafTools require legacy Python (Docker fallback stronglyrecommended)

## Future Development

- Complete alignment-free rule integration
- Improve checkpoint/resume functionality
- Add dynamic resource allocation
- Optimize memory usage for large genomes
- Parallelize annotation steps further

For detailed script documentation, see README files in `scripts/step_*/` subdirectories.