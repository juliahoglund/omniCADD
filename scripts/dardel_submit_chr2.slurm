#!/bin/bash
#SBATCH -A naiss2024-22-1246
#SBATCH -p main
#SBATCH -t 04:00:00
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -J omniCADD_chr2_snpeff

set -euo pipefail

# Load modules as needed on Dardel
module load python3 || true
module load conda || true

# Optional: point conda envs to project space
# export CONDA_ENVS_PATH=/cfs/klemming/projects/snic/snic2022-22-894/condaenvs

# Run Snakemake orchestrator inside this Slurm job and submit worker jobs via the Slurm profile
snakemake -s run_chr2_test.smk \
  --configfile config/profile_limited_data.yaml \
  --profile config/slurm/profile \
  --use-conda \
  --jobs 50 \
  --printshellcmds \
  target_chr2
