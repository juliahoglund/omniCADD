#!/usr/bin/env bash
# Minimal HPC runner for SNPEff chr2 test
# Usage: bash scripts/run_hpc_chr2.sh

set -euo pipefail

# Optional: load modules (example placeholders)
# module load bioinfo-tools
# module load conda

# Run with conda-managed envs
snakemake -s run_chr2_test.smk \
  --configfile config/profile_limited_data.yaml \
  --use-conda \
  --cores 4 \
  --printshellcmds \
  target_chr2
