#!/bin/bash
#SBATCH -A naiss2025-5-78
#SBATCH -p main
#SBATCH -t 06:00:00
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -J omniCADD_augustus_snpeff_tests
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.hoglund@su.se

set -euo pipefail

module load python3 || true
module load snakemake/7.25.0 || true
conda config --set channel_priority strict

# Fallback: use user Miniforge if conda module is unavailable
if [ -d "$HOME/miniforge3" ]; then
  export PATH="$HOME/miniforge3/bin:$PATH"
  # shellcheck disable=SC1091
  source "$HOME/miniforge3/etc/profile.d/conda.sh" 2>/dev/null || true
fi

# Optional: conda envs path in project space
# export CONDA_ENVS_PATH=/proj/naiss2025-5-78/condaenvs

# Preflight against Augustus overlay
python3 scripts/preflight_check.py --profile config/test_overlay_augustus.yaml || { echo "Preflight failed"; exit 1; }

# Run multi-chromosome test wrapper (chr2 + chr18)
snakemake -s run_multi_test.smk \
  --configfile config/test_overlay_augustus.yaml \
  --profile config/slurm/profile \
  --use-conda \
  --conda-frontend mamba \
  --jobs 100 \
  --printshellcmds \
  target_multi
