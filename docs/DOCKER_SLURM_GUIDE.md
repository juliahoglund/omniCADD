# Docker and SLURM Configuration Guide

## Changes Made

### 1. Docker Containers Activated
- **Uncommented** all mafTools container references in `workflow/rules/1_extract_ancestor.smk`
- **Fixed** Dockerfile best practices:
  - SIFT4G: Uses specific Ubuntu version (20.04) and proper cleanup
  - renv: Fixed hardcoded path reference
- **Created** build script: `build_containers.sh`

### 2. SLURM Configuration Updated

#### Profile Configuration (`config/slurm/profile/config.yaml`):
- **Enabled** Singularity support (`use-singularity: True`)
- **Added** proper Singularity arguments for container isolation
- **Configured** job control parameters for cluster stability

#### Cluster Resources (`config/slurm/cluster.yaml`):
- **High-resource rules**: mafTools operations get 6-8 CPUs, 8GB RAM
- **Memory-intensive rules**: GERP annotation gets 8 CPUs, 16GB RAM
- **Model training**: 16 CPUs, 64GB RAM, 24 hours
- **Scoring**: 8 CPUs, 32GB RAM, 12 hours

## How to Use

### Building Docker Containers

1. **Build all containers:**
   ```bash
   ./build_containers.sh
   ```

2. **Build specific container:**
   ```bash
   ./build_containers.sh mafTools
   ```

3. **Push to Docker Hub** (after building):
   ```bash
   ./build_containers.sh push
   ```

### Running with SLURM

1. **Local testing:**
   ```bash
   snakemake --use-conda --use-singularity --cores 4
   ```

2. **SLURM cluster:**
   ```bash
   snakemake --profile config/slurm/profile --cluster-config config/slurm/cluster.yaml
   ```

3. **Specific SLURM system (e.g., Rackham):**
   ```bash
   snakemake --use-conda --use-singularity --profile config/slurm/profile \
             --cluster-config config/slurm/cluster.yaml \
             --conda-prefix /path/to/your/conda/envs
   ```

## Configuration Notes

### Container Strategy
- **mafTools**: Docker container (Python 2.7 dependency)
- **GERP**: Singularity from Biocontainers
- **Other tools**: Conda environments (faster, easier to manage)

### Resource Allocation
- **Default**: 1 CPU, 5 hours, minimal memory
- **mafTools rules**: 2-6 CPUs, 8GB RAM, 6-8 hours
- **Annotation**: 8 CPUs, 16GB RAM, 12 hours
- **Training**: 16 CPUs, 64GB RAM, 24 hours

### Troubleshooting

1. **Container not found**: Build containers first with `./build_containers.sh`
2. **Singularity issues**: Check cluster has Singularity installed
3. **Memory errors**: Increase memory limits in `cluster.yaml`
4. **Time limits**: Adjust time in `cluster.yaml` based on your data size

### Customization

1. **Change Docker registry**: Edit `REGISTRY_PREFIX` in `build_containers.sh`
2. **Adjust resources**: Modify `config/slurm/cluster.yaml`
3. **Add new rules**: Follow the pattern in `cluster.yaml`
4. **Local conda prefix**: Set in your SLURM profile or command line
