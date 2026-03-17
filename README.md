# omniCADD
#### A CADD scoring system to assess variant deleteriousness in non-model organisms
----

[![Documentation Status](https://readthedocs.org/projects/omnicadd/badge/?version=latest)](https://omnicadd.readthedocs.io/en/latest/?badge=latest)

## Overview

omniCADD is a **Snakemake-based bioinformatics pipeline** that computes CADD (Combined Annotation Dependent Depletion) scores—a measure of variant deleteriousness—for **non-model organisms**. This repository contains a reconstructed workflow for computing CADD scores, as first described in [Kircher et al 2014](https://www.nature.com/articles/ng.2892). 

The pipeline combines:
- **Evolutionary information** from multi-species alignments
- **Functional annotations** (VEP/SNPEff)
- **Conservation metrics** (PHAST, GERP)
- **Machine learning** (logistic regression)

To generate **PHRED-scaled deleteriousness scores** for every possible variant in a genome.

The code is an update, modification and extension based on the CADD models for [mouse](https://doi.org/10.1186/s12859-018-2337-5) (Groß et al 2018), [pig](https://doi.org/10.1186/s12711-020-0528-9) (Groß et al 2020a), and [chicken](https://doi.org/10.1371/journal.pgen.1009027) (Groß et al 2020b).

## Workflow Overview

omniCADD processes genomic data through 8 main steps:

1. **Extract Ancestral Sequences** - Reconstruct ancestral sequences from multi-species alignments
2. **Derive Variants** - Identify derived variants from population VCF data
3. **Simulate Variants** - Generate variants for model training
4. **Summary Report** - Create HTML visualization of summary statistics of step 1, 2 and 3 * 
5. **Annotate Variants** - Apply functional (VEP) and conservation (PHAST, GERP) annotations
6. **Combine Annotations** - Merge all annotations into feature matrices
7. **Train & Test Model** - Build logistic regression classifier with cross-validation
8. **Score Variants** - Generate genome-wide PHRED-scaled CADD scores

\* Additional types of annotations are under development 

**Training Data**:
- **Derived variants** (Step 2) from population data → labeled as **benign** (class 0)
  - Variants observed at appreciable frequencies have passed through natural selection
- **Simulated variants** (Step 3) randomly generated → labeled as **deleterious** (class 1)
  - Random mutations not filtered by selection, enriched for harmful variants

See the [workflow documentation](workflow/README.md) for detailed information about each step.

## Pipeline Modes

### Main Pipeline - Dense Data
As of now, the main branch contains the workflow for running it using already deposited data, such as EPO alignments from Ensembl, deposited reference genomes, and Ensembl VEP. It also assumes that the user has a population-wide VCF file with variants to be used in the section *derive variants*. 

### Alignment-Free Pipeline - Scarce Data
The *alignmentfree* branch (still under development) assumes the opposite, namely that:

- The species of interest is not part of a publicly available alignment
- The reference species is not deposited in Ensembl
- The species of interest does not need to have gene annotation (gff/gff3)

The plan is to merge and make the pipeline fully modifiable in the future. 

#### :exclamation: **NOTE** :exclamation: 
The alignment-free scripts are still under construction and have not yet been fully integrated into the workflow. 

## Requirements

### Software Dependencies
- **Snakemake** ≥ 7.0
- **Conda/Mamba** for environment management
- **Docker** or **Singularity/Apptainer** (optional, for containerized execution)*

\* the plan is to make it fully containerised in the future

### System Requirements
- **RAM**: 16+ GB recommended (64+ GB for whole genomes)
- **Storage**: 500+ GB for alignments, annotations, and intermediate files
- **CPU**: Multi-core processor recommended for parallelization

### Conda Environments
The pipeline uses multiple specialized conda environments (see [workflow/envs/](workflow/envs/)):
- `ancestor.yml` - Ancestral sequence processing
- `annotation.yml` - VEP, PHAST, GERP conservation tools
- `model.yml` - Machine learning (scikit-learn)
- `score.yml` - Variant generation and scoring
- `simulation.yml` - Variant simulation
- `common.yml` - Utility tools (bgzip, tabix, bcftools)
- `r_stats.yml` - R environment for reporting

## Input Data Requirements

### Required Inputs
The pipeline expects the following data in the `resources/` directory:

| Data Type | Location | Description | Example Source |
|-----------|----------|-------------|----------------|
| **Reference Genome** | `resources/genome/` | FASTA format, gzipped | Ensembl FTP |
| **Multi-species Alignment** | `resources/alignment/` | MAF format (EPO alignments) | Ensembl Compara |
| **Population VCF** | `resources/pop-level-VCF/` | VCF files with population variants | Your data |
| **Phylogenetic Tree** | `resources/` | Newick format (.nwk) | Ensembl Compara |

### Optional Inputs
| Data Type | Location | Purpose |
|-----------|----------|---------|
| **Gene Annotation** | `resources/` | GFF3/GTF format for visualization |


See [resources/README.md](resources/README.md) for detailed data preparation instructions.

## :link: Quick Start 

### 1. Download Required Data
To work correctly, the pipeline assumes that critical files are in subfolders in `resources/`. Download or transfer them before running:

#### Reference Genome
```bash
cd resources/genome
# Example: Download pig reference genome. (This is what has been used to build the pipeline)
wget https://ftp.ensembl.org/pub/current_fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
```

#### Multi-species Alignment
```bash
cd ../alignment
# Example: Download 44-way mammal EPO alignment (This is what has been used to build the pipeline)
wget -r -nd --no-parent -e robots=off -A '44_mammals.epo.*.maf.gz' \
  https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/44_mammals.epo/
```

#### Population VCF
```bash
cd ../pop-level-VCF
# Transfer your population VCF files here
# Files should contain variants with allele frequencies
```

#### Phylogenetic Tree
```bash
cd ../
# Example: Download 44 mammal tree (This is what has been used to build the pipeline)
wget https://ftp.ensembl.org/pub/current_compara/species_trees/44_eutherian_mammals_EPO_default.nh
mv 44_eutherian_mammals_EPO_default.nh tree_43_mammals.nwk
```

#### Optional: Gene Annotations
```bash
# GFF3 and GTF files (for visualization in summary report)
# These are optional and only affect cosmetic aspects of the report
```

### 2. Configure the Pipeline
Edit `config/config.yaml` to specify:
- Species name and chromosomes
- Ancestor species and comparison species
- Annotation methods (VEP, SNPEff, or both)
- Model parameters

See [config/README.md](config/README.md) for configuration details.

### 3. Run the Pipeline
```bash
# Dry run to check workflow
snakemake -n

# Local execution with 8 cores
snakemake --use-conda --cores 8

# SLURM cluster execution
snakemake --profile config/slurm/profile/
```

## Output

The pipeline generates PHRED-scaled CADD scores for all possible variants:

**Main Output**: `results/cadd_scores/chr*.tsv.gz`
- Tab-separated, bgzip-compressed, tabix-indexed
- Format: `CHROM  POS  REF  ALT  CADD_SCORE`
- PHRED scores: higher = more deleterious

**Intermediate Outputs**:
- `results/ancestral_seq/` - Reconstructed ancestor sequences
- `results/derived_variants/` - Real variants from population (benign training set)
- `results/simulated_variants/` - Simulated variants (deleterious-enriched training set)
- `results/annotations/` - VEP, PHAST, GERP annotations
- `results/models/` - Trained logistic regression classifier
- `output/*.html` - Summary reports

## Repository Structure

```
omniCADD/
├── snakefile             # Main Snakemake workflow entry point
├── config/               # Configuration files
│   ├── config.yaml       # Main pipeline configuration
│   ├── *.tsv             # Annotation combination configs
│   └── slurm/            # SLURM cluster profiles
├── workflow/             # Workflow components
│   ├── rules/            # Snakemake rules (8 workflow steps)
│   ├── scripts/          # Python/Perl/R scripts for each step
│   ├── envs/             # Conda environment definitions
│   └── docker/           # Dockerfiles
├── resources/            # Input data directory
│   ├── genome/           # Reference genome (FASTA)
│   ├── alignment/        # Multi-species alignments (MAF)
│   └── pop-level-VCF/    # Population variant data
├── results/              # Output directory (auto-generated)
└── output/               # Summary reports
```

See individual README files in each directory for more details.

## Current Limitations

1. **Population VCF Required**: Variant derivation requires high-frequency variants (>90% AF)
2. **EPO Alignments**: Main branch assumes species is in Ensembl EPO alignments
3. **Python 2.7 Dependency**: Some mafTools require legacy Python (wrapped in Docker)
4. **Alignment-Free Branch**: Still under development, not fully integrated
5. **Memory Requirements**: Large alignments may require 64+ GB RAM
6. **Validation Data**: No built-in support for external validation sets

## Future Development

### High Priority
- [ ] **Complete alignment-free branch** for poorly-annotated species
- [ ] **Variant derivation without population VCF** using alternative methods
- [ ] **Merge standard and alignment-free pipelines** into unified modifiable workflow
- [ ] **Add validation set support** for model evaluation with known pathogenic variants

### Pipeline Improvements
- [ ] Clean up and refactor scripts for better maintainability
- [ ] Update gene annotation handling for better visualization
- [ ] Make cluster profiles more production-ready
- [ ] Add checkpoint/resume functionality for long-running jobs
- [ ] Implement dynamic resource allocation based on data size

### Feature Enhancements
- [ ] Support for indel scoring (currently basic)
- [ ] Integration with ab initio gene prediction tools (Augustus)
- [ ] Multi-species preset configurations
- [ ] Improved SIFT score handling when available

In a galaxy far far away, maybe some day: 
- [ ] Web interface for score queries
- [ ] Downstream analysis tools (variant prioritization, cohort analysis)

### Documentation
- [ ] Best practices guide for different species
- [ ] Performance benchmarks

## Information

See [References](#references) section for studies of which this pipeline is based and adapted from.

### Support & Contributing

- **Issues**: Report bugs or request features via GitHub Issues
- **Documentation**: Full documentation at https://omnicadd.readthedocs.io/
- **Contributing**: Pull requests welcome! See CONTRIBUTING.md (coming soon)

## :ballot_box_with_check: TO-DO

### Critical
- Clean up scripts and workflow
- Remove temporary debugging code
- Implement running without population-level VCF files
- Finish alignment-free branch integration

### Enhancements
- Make cluster profile more robust
- Alter visualization to work with in-pipeline ab initio gene predictions
- Add comprehensive testing

## References
Kircher, M., Witten, D., Jain, P. et al. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet 46, 310–315 (2014). doi.org/10.1038/ng.2892

Groß, C., de Ridder, D. & Reinders, M. Predicting variant deleteriousness in non-human species: applying the CADD approach in mouse. BMC Bioinformatics 19, 373 (2018). doi.org/10.1186/s12859-018-2337-5

Groß, C., Derks, M., Megens, HJ. et al. pCADD: SNV prioritisation in Sus scrofa. Genet Sel Evol 52, 4 (2020). doi.org/10.1186/s12711-020-0528-9

Groß C, Bortoluzzi C, de Ridder D, Megens HJ, Groenen MAM, et al. (2020) Prioritizing sequence variants in conserved non-coding elements in the chicken genome using chCADD. PLOS Genetics 16(9): e1009027. doi.org/10.1371/journal.pgen.1009027

McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F.
The Ensembl Variant Effect Predictor.
Genome Biology Jun 6;17(1):122. (2016)
doi:10.1186/s13059-016-0974-4 