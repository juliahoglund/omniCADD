# omniCADD Setup Guide for PDC Dardel

This guide provides specific setup instructions for running omniCADD on PDC Dardel.

## Prerequisites on Dardel

### 1. Load Required Modules
```bash
module load PDC bioinfo-tools apptainer
```

### 2. Set Up Storage Paths

```bash
# Set cache directories (add to ~/.bashrc)
export PROJ=/cfs/klemming/projects/supr/YOUR_PROJECT
export XDG_CACHE_HOME=$PROJ/omniCADD/.cache
export APPTAINER_CACHEDIR=$PROJ/apptainer-cache

# Create directories
mkdir -p $XDG_CACHE_HOME $APPTAINER_CACHEDIR
```

## Setup Steps

### 1. Clone and Setup Repository
```bash
cd /cfs/klemming/projects/supr/YOUR_PROJECT/
git clone <your-omniCADD-repo>
cd omniCADD
mkdir -p slurm   # directory for SLURM job logs
```

### 2. Create Conda Environment
```bash
# Keep conda environments in project storage
export CONDA_ENVS_PATH=/cfs/klemming/projects/supr/YOUR_PROJECT/conda-envs/

conda env create -f workflow/envs/common.yml -p $CONDA_ENVS_PATH/omnicadd
conda activate $CONDA_ENVS_PATH/omnicadd

# Install the SLURM executor plugin (required for Snakemake 8+)
pip install snakemake-executor-plugin-slurm
```

### 3. Configure the Dardel Profile

Edit `config/slurm/dardel/config.yaml` and fill in the two placeholders:

```yaml
# Your PDC allocation (shown in `projinfo` or the PDC portal, e.g. naiss2024-X-NNN)
slurm_account: "YOUR_DARDEL_ACCOUNT"

# Absolute path to your project storage for apptainer image cache
apptainer-prefix: "/cfs/klemming/projects/supr/YOUR_PROJECT/apptainer_cache"
```

### 4. Update Main Config

Edit `config/config.yaml` to point to your data on project storage:
```yaml
# Use project storage for large files
generate_variants:
  reference_genome_wildcard: "/cfs/klemming/projects/supr/YOUR_PROJECT/..."
```

## Running on Dardel

### 1. Start a Screen Session (keep running after logout)
```bash
screen -S omnicadd

# Or with tmux:
tmux new -s omnicadd
```

### 2. Activate Environment and Load Modules
```bash
module load PDC bioinfo-tools apptainer
export CONDA_ENVS_PATH=/cfs/klemming/projects/supr/YOUR_PROJECT/conda-envs/
conda activate $CONDA_ENVS_PATH/omnicadd
cd /cfs/klemming/projects/supr/YOUR_PROJECT/omniCADD
```

### 3. Test Setup (Dry Run)
```bash
snakemake --profile config/slurm/dardel -n
```

### 4. Run Pipeline
```bash
# Full run — Snakemake stays on the login node and submits each rule as a separate SLURM job
snakemake --profile config/slurm/dardel

# Re-run incomplete jobs and keep going on errors:
snakemake --profile config/slurm/dardel --rerun-incomplete --keep-going
```

SLURM job logs are written to `slurm/slurm_<jobid>.out` / `.err`.  
To monitor running jobs: `squeue -u $USER`

## Dardel-Specific Optimizations

### 1. Container Strategy
- **Apptainer** (formerly Singularity) is used instead of Docker
- Containers are automatically pulled/converted by Snakemake
- Cache stored in project storage to avoid home directory limits

### 2. Environment Strategy
- Conda environments are stored in project storage to avoid home directory quota limits
- Each Snakemake rule uses its own conda environment defined in `workflow/envs/`

### 4. Storage Strategy
- **Project storage** for all large files and caches
- **Home directory** only for configuration files
- **Scratch space** ($PDC_TMP) for temporary files

### Useful Dardel Commands
```bash
# Check project storage usage
projinfo <project-id>

# Check job status
squeue -u $USER

# Check detailed job info
scontrol show job <job-id>
```