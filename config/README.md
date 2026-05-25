# Configuration Directory

This directory contains all configuration files for the omniCADD pipeline.

## Directory Structure

```
config/
├── config.yaml                    # Main pipeline configuration
├── annot_combinations_config.tsv  # Annotation combination settings  
├── annot_processing_config.tsv    # Annotation processing parameters
├── slurm/                         # SLURM cluster profiles
│   ├── profile/                   # Snakemake profiles
│   ├── rackham/                   # UPPMAX Rackham configs
│   └── dardel/                    # PDC Dardel configs
└── README.md                      # This file
```

## Main Configuration: `config.yaml`

The primary configuration file controls all aspects of the pipeline.

### Species Configuration
Configuration examples stated below is how the pipeline was build and tested. 

```yaml
species_name: "sus_scrofa"  # Reference species
comparison: "cow"           # For ancestral reconstruction
```

**Important**: Species names must match those in alignment and phylogenetic tree exactly.

### Chromosome Configuration

```yaml
chromosomes:
  autosomes: [1, 2, 3, ..., 18]
  karyotype: [1, 2, 3, ..., 18, "X", "Y"]
  train: [1, 2, 3, 4, 5]        # Chromosomes for model training
  score: [1, 2, 3, 4, 5]        # Chromosomes to score
```

**Settings**:
- `autosomes`: All autosomal chromosomes
- `karyotype`: Complete karyotype (including sex chromosomes)
- `train`: Subset for model training
- `score`: Chromosomes to generate final CADD scores for

It is recommended to create genome wide-scores to get genome-wide phred scaling.
However, a smaller subset can be used for testing workflow and param settings before running the full analysis.

### Ancestral Sequence Mode

```yaml
ancestor:
  mode: "alignment"  # Options: alignment | reconstruct
```

**Modes**:
- **`alignment`** (Recommended): Extract from EPO multi-species alignments
- **`reconstruct`** (Under development): Align outgroup genome to reference

### Annotation Configuration

```yaml
annotation:
  variant_annotator: "vep"  # Options: "vep", "snpeff", "both"
  vep_version: 111          # Ensembl VEP cache version
  vep_auto_install: true    # Auto-download VEP cache
  
  conservation:
    phast: true             # Run PHAST conservation
    gerp: true              # Run GERP conservation
```

**Variant Annotators**:
- `vep`: Ensembl Variant Effect Predictor (recommended when available)
- `snpeff`: SNPEff annotation
- `both`: Run both VEP and SNPEff (slower, more comprehensive)

### Model Parameters

```yaml
model:
  algorithm: "logistic_regression"
  k_folds: 5  # Cross-validation folds
  
  grid_search:
    C_values: [0.001, 0.01, 0.1, 1, 10]    # Regularization
    max_iter: [1000, 2000, 5000]           # Max iterations
```

Pipeline tests all combinations of C and max_iter to find optimal hyperparameters.

### Variant Configuration

```yaml
derived_variants:
  af_threshold: 0.90          # Minimum allele frequency (90%)
  quality_filter: "QUAL>30"   # VCF quality filter

simulated_variants:
  overestimation_factor: 5    # Mutation rate scaling
  simulate_indels: false      # Include INDELs (slower, scoring not yet implemented)
```

---

## SLURM Cluster Configuration

### Cluster Resource Configuration

Each workflow step has custom resource requirements in `cluster.yaml`:

```yaml
annotate_vep:
  account: ""
  partition: "core"
  time: "12:00:00"
  cpus: 4
  mem: "32GB"

train_model:
  account: ""
  partition: "core"
  time: "24:00:00"
  cpus: 8
  mem: "64GB"
```

**Customize** these for your cluster allocation and requirements.

---

## Container Usage

### Run with Singularity/Apptainer (HPC Clusters)

For tools difficult to install (e.g., mafTools, R environment):

```bash
# Snakemake rule example:
rule extract_ancestor:
    container:
        "docker://juliahoglund/maftools:latest"
    # ... rest of rule

# Run pipeline with Singularity
snakemake --use-singularity --use-conda
```

### Available Containers

| Container | Purpose | Image |
|-----------|---------|-------|
| mafTools | Python 2.7 MAF tools | `juliahoglund/maftools:latest` |
| GERP | GERP conservation | `juliahoglund/gerp:latest` |
| SIFT4G | SIFT4G scoring | `juliahoglund/sift4g:latest` |
| R environment | Complete R env with all packages | `juliahoglund/report-renv:latest` |

## Annotation Configuration Files

### `annot_combinations_config.tsv`

Defines which annotation combinations to test:

```tsv
Combination	VEP	PHAST	GERP
full	TRUE	TRUE	TRUE
vep_only	TRUE	FALSE	FALSE
conservation_only	FALSE	TRUE	TRUE
```

Use this to evaluate feature importance by testing subsets of annotations.

### `annot_processing_config.tsv`

Parameters for annotation processing:

```tsv
Parameter	Value
chunk_size	10000
min_alignment_depth	3
max_gap_fraction	0.5
```

---

## Troubleshooting

### Configuration Validation

```bash
# Check YAML syntax
python -c "import yaml; yaml.safe_load(open('config/config.yaml'))"

# Dry run to validate configuration
snakemake -n
```

### Common Issues

**Issue**: Species name mismatch  
**Solution**: Ensure `species_name` matches alignment and tree exactly (case-sensitive)

**Issue**: VEP cache download fails  
**Solution**: Set `vep_auto_install: false` and install manually

**Issue**: Cluster jobs fail immediately  
**Solution**: Check account, partition, and resource limits in cluster config

**Issue**: Out of memory  
**Solution**: Increase memory in cluster config or reduce chunk_size

---

## Configuration Templates

Example configurations in repository:

- `config/config.yaml` - Default pig configuration
- `config/slurm/profile/config.yaml` - UPPMAX profile
- `config/slurm/profile/config_dardel.yaml` - Dardel profile

Copy and modify for your species and cluster environment.

## Further Reading

- [Snakemake Profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
- [SLURM Documentation](https://slurm.schedmd.com/sbatch.html)
- [Singularity User Guide](https://sylabs.io/guides/latest/user-guide/)
- [Main README](../README.md) for overall setup instructions
