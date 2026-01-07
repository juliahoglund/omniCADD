#!/bin/bash
#SBATCH -A naiss2025-5-78
#SBATCH -p main
#SBATCH -t 04:00:00
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -J omniCADD_tests_overlay

set -euo pipefail

module load python3 || true
module load conda || true
module load singularity || true

# Optional: conda envs path
# export CONDA_ENVS_PATH=/cfs/klemming/projects/snic/snic2022-22-894/condaenvs

# Preflight against overlay
python3 scripts/preflight_check.py --profile config/test_overlay.yaml || { echo "Preflight failed"; exit 1; }

# Run with layered configs: base + overlay
snakemake -s run_multi_test.smk \
  --configfile config/test_overlay.yaml \
  --profile config/slurm/profile \
  --use-conda \
  --use-singularity \
  --singularity-args "-B $PWD" \
  --jobs 50 \
  --printshellcmds \
  target_multi
