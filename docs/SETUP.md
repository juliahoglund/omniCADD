# omniCADD Setup Guide

This guide helps you set up the omniCADD pipeline for different execution environments.

## Prerequisites

- Snakemake (≥8.0.0 recommended; ≥7.21.0 via conda also works)
- Conda/Mamba
- Docker (for local containerized execution)
- Singularity (for cluster containerized execution)

## Quick Start

### Local Execution

```bash
# Create conda environments
snakemake --use-conda --conda-create-envs-only

# Dry run
snakemake -n

# Full run
snakemake --cores 8 --use-conda
```

### Cluster Execution (SLURM)

**For UPPMAX (Rackham/Snowy/Bianca):**
1. **Update conda prefix** in `config/slurm/profile/config.yaml`:
   ```yaml
   conda-prefix: "/path/to/your/conda/envs"
   ```

2. **Run with cluster profile**:
   ```bash
   # Create cluster-specific environments
   snakemake --use-conda --conda-create-envs-only --profile config/slurm

   # Submit jobs
   snakemake --profile config/slurm
   ```

**For PDC Dardel:**
1. **Use Dardel-specific configuration**:
   ```bash
   # Copy Dardel configuration
   cp config/slurm/profile/config_dardel.yaml config/slurm/profile/config.yaml
   cp config/slurm/cluster_dardel.yaml config/slurm/cluster.yaml
   ```

2. **Load required modules and run**:
   ```bash
   module load PDC bioinfo-tools apptainer tmux
   snakemake --profile config/slurm --directory /path/to/omniCADD/
   ```

   See `DARDEL_SETUP.md` for detailed Dardel-specific instructions.

## Container Setup

### For Docker (Local)
```bash
# Build all containers
bash workflow/build_containers.sh
```

### For Singularity (Cluster)
The pipeline will automatically pull/convert containers when needed.

## Environment Switching

The pipeline automatically detects your environment:
- **Local**: Uses generic conda environments from `workflow/envs/`
- **Cluster**: Uses cluster-specific environments from `config/slurm/rackham/`

## Configuration

### Main Config (`config/config.yaml`)
- Set input/output paths
- Configure alignment and annotation parameters
- Adjust species and chromosome settings

### Cluster Config (`config/slurm/cluster.yaml`)
- Adjust memory/time limits per rule
- Configure partition preferences

## Troubleshooting

### Conda Environment Issues
- Ensure conda prefix is correctly set in cluster profile
- Create environments manually if needed:
  ```bash
  snakemake --use-conda --conda-create-envs-only
  ```

### Container Issues
- For Docker: Check daemon is running
- For Singularity: Ensure proper bind paths in cluster config

### Memory/Time Limits
- Adjust resource requirements in `config/slurm/cluster.yaml`
- Monitor job logs for specific failures

## File Structure

```
omniCADD/
├── config/                    # Configuration files
│   ├── config.yaml           # Main pipeline config
│   └── slurm/                # Cluster-specific configs
├── workflow/                 # Snakemake workflow
│   ├── rules/                # Rule definitions
│   ├── envs/                 # Generic conda environments
│   ├── scripts/              # Helper scripts
│   └── docker/               # Container definitions
└── results/                  # Output directory (created during run)
```

For detailed information about specific components, see:
- `workflow/README.md` - Workflow details
- `config/README.md` - Configuration options
