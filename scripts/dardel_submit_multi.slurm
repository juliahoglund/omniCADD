#!/bin/bash
#SBATCH -A naiss2025-5-78
#SBATCH -p main
#SBATCH -t 06:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -J omniCADD_multi_snpeff

set -euo pipefail

# Load modules as needed on Dardel
module load python3 || true
module load conda || true
module load singularity || true

# Use dedicated cache in project workspace to avoid HOME quota
export SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-"$PWD/.singularity/cache"}
export SINGULARITY_TMPDIR=${SINGULARITY_TMPDIR:-"$PWD/.singularity/tmp"}
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" || true

# Keep conda package cache out of HOME to reduce quota pressure
export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS:-"$PWD/.conda/pkgs"}
mkdir -p "$CONDA_PKGS_DIRS" || true

# Optional: point conda envs to project space
# export CONDA_ENVS_PATH=/cfs/klemming/projects/snic/snic2022-22-894/condaenvs

# Preflight validation (optional; comment out if not desired)
python3 scripts/preflight_check.py --profile config/profile_chr2_chr18.yaml || { echo "Preflight failed"; exit 1; }

# Run Snakemake orchestrator and submit worker jobs via the Slurm profile
snakemake -s run_multi_test.smk \
  --configfile config/profile_chr2_chr18.yaml \
  --profile config/slurm/profile \
  --use-conda \
  --use-singularity \
  --singularity-args "-B $PWD" \
  --latency-wait 60 \
  --conda-prefix "$PWD/.snakemake/conda" \
  --rerun-incomplete \
  --keep-going \
  --jobs 100 \
  --printshellcmds \
  target_multi
