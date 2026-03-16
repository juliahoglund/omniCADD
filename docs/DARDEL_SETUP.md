# omniCADD Setup Guide for PDC Dardel

This guide provides specific setup instructions for running omniCADD on PDC Dardel.

## Prerequisites on Dardel

### 1. Load Required Modules
```bash
module load PDC bioinfo-tools apptainer tmux
```

### 2. Set Up Storage Paths

**Important**: Dardel has limited home directory space. Use your project storage:

```bash
# Set cache directories (add to ~/.bashrc)
export XDG_CACHE_HOME=/cfs/klemming/projects/supr/sllstore.../omniCADD/.cache
export APPTAINER_CACHEDIR=/cfs/klemming/projects/supr/sllstore.../apptainer-cache

# Create directories
mkdir -p $XDG_CACHE_HOME
mkdir -p $APPTAINER_CACHEDIR
```

## Setup Steps

### 1. Clone and Setup Repository
```bash
cd /cfs/klemming/projects/supr/sllstore.../
git clone <your-omniCADD-repo>
cd omniCADD
chmod 755 snakefile
```

### 2. Create Conda Environment
```bash
# Set conda environment path to project storage
export CONDA_ENVS_PATH=/cfs/klemming/projects/supr/sllstore.../conda-envs/

# Create environment
conda env create -f workflow/envs/omnicadd.yml -p $CONDA_ENVS_PATH/omnicadd
```

### 3. Configure for Dardel

Copy and edit the Dardel-specific configuration:
```bash
# Copy Dardel SLURM profile
cp config/slurm/profile/config_dardel.yaml config/slurm/profile/config.yaml

# Copy Dardel cluster configuration  
cp config/slurm/cluster_dardel.yaml config/slurm/cluster.yaml
```

**Edit paths in both files:**
- Update `conda-prefix` with your actual project path
- Update `singularity-prefix` with your actual project path
- Add your PDC project account in `cluster.yaml`

### 4. Update Main Config

Edit `config/config.yaml`:
```yaml
# Use project storage for large files
input_dir: "/cfs/klemming/projects/supr/sllstore.../data"
output_dir: "/cfs/klemming/projects/supr/sllstore.../results"
```

## Running on Dardel

### 1. Start Interactive Session
```bash
# Start tmux session
tmux new -s omnicadd

# Activate conda environment
export CONDA_ENVS_PATH=/cfs/klemming/projects/supr/sllstore.../conda-envs/
conda activate omnicadd
```

### 2. Test Setup (Dry Run)
```bash
snakemake --profile config/slurm -n --directory /cfs/klemming/projects/supr/sllstore.../omniCADD/
```

### 3. Run Pipeline
```bash
# Full run
snakemake --profile config/slurm --directory /cfs/klemming/projects/supr/sllstore.../omniCADD/

# With useful flags
snakemake --profile config/slurm --ri -k --directory /cfs/klemming/projects/supr/sllstore.../omniCADD/
```

## Dardel-Specific Optimizations

### 1. Container Strategy
- **Apptainer** (formerly Singularity) is used instead of Docker
- Containers are automatically pulled/converted by Snakemake
- Cache stored in project storage to avoid home directory limits

### 2. Resource Allocation
- Dardel nodes have 128 cores and 256GB RAM
- Configuration optimized for these specifications
- Uses `main` partition for general jobs

### 3. Environment Strategy
- **Dardel-specific** conda environments in `config/slurm/dardel/`
- **Fallback** to generic environments if Dardel-specific not found
- **Intel-optimized** packages where available

### 4. Storage Strategy
- **Project storage** for all large files and caches
- **Home directory** only for configuration files
- **Scratch space** ($PDC_TMP) for temporary files

## Troubleshooting

### Common Issues

1. **Home directory full**: Move all caches to project storage
2. **Module conflicts**: Use exact module versions shown above  
3. **Container pulls fail**: Check apptainer cache permissions
4. **Jobs killed**: Increase memory in `cluster.yaml`

### Performance Tips

1. **Use Intel-optimized packages** where available
2. **Optimize I/O** by using fast project storage
3. **Monitor resource usage** and adjust cluster config
4. **Use scratch space** for temporary files when possible

### Useful Dardel Commands
```bash
# Check project storage usage
projinfo <project-id>

# Check job status
squeue -u $USER

# Check detailed job info
scontrol show job <job-id>
```

For more details about Dardel, see: https://www.pdc.kth.se/hpc-services/computing-systems/dardel
