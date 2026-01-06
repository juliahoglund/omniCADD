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

# Reset modules to avoid pulling old Rackham paths
module purge || true
ml PDC || true
ml snakemake || true

# Load container runtime (force apptainer; unload singularity to avoid Snakemake picking it)
ml apptainer || true
module unload singularity >/dev/null 2>&1 || true

# Verify container runtime availability
if ! command -v apptainer >/dev/null 2>&1; then
  echo "ERROR: apptainer/singularity not available on this node. Load the appropriate module (e.g., 'module load apptainer')." >&2
  exit 1
fi

# Select runtime binary (apptainer only)
CONTAINER_RUNTIME="apptainer"

# Use dedicated cache in HOME (or project scratch if preferred)
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-"$HOME/.apptainer/cache"}
export SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-"$APPTAINER_CACHEDIR"}
mkdir -p "$APPTAINER_CACHEDIR" || true


# Prefer user Miniforge/conda; configure strict channel priority
if [ -d "$HOME/miniforge3" ]; then
  export PATH="$HOME/miniforge3/bin:$PATH"
  # shellcheck disable=SC1091
  source "$HOME/miniforge3/etc/profile.d/conda.sh" 2>/dev/null || true
fi
conda config --set channel_priority strict || true

# Disable conda plugins to avoid environment issues
export CONDA_NO_PLUGINS=true
# Force classic solver to avoid libmamba plugin requirement
export CONDA_SOLVER=classic

# Use classic conda frontend for reliability on Dardel
FRONTEND="conda"

# Optional: conda envs path in project space
# export CONDA_ENVS_PATH=/proj/naiss2025-5-78/condaenvs

# Preflight against Augustus overlay
python3 scripts/preflight_check.py --profile config/test_overlay_augustus.yaml || { echo "Preflight failed"; exit 1; }

rm -rf .snakemake/conda/287305f04c23dece5b5947f7c78de31d_

# Helper: detect Augustus config dir inside container
detect_cfg() {
  local img="$1"
  for p in \
    "/opt/conda/envs/augustus/share/augustus/config" \
    "/opt/conda/share/augustus/config" \
    "/opt/conda/config/augustus" \
    "/usr/share/augustus/config" \
    "/usr/local/share/augustus/config"; do
    if $CONTAINER_RUNTIME exec "$img" bash -lc "[ -d '$p' ]"; then
      echo "$p"
      return 0
    fi
  done
  # Try environment variable inside container
  local envp
  envp=$($CONTAINER_RUNTIME exec "$img" bash -lc 'echo ${AUGUSTUS_CONFIG_PATH:-}' 2>/dev/null || true)
  if [ -n "$envp" ] && $CONTAINER_RUNTIME exec "$img" bash -lc "[ -d '$envp' ]"; then
    echo "$envp"
    return 0
  fi
  return 1
}

# Quick sanity: extract Augustus image from overlay and verify only if it's a local SIF
AUG_IMG=$(awk '/^containers:/{f=1} f && /augustus:/ {print; exit}' config/test_overlay_augustus.yaml | sed -E 's/.*"(.*)".*/\1/')
if [ -n "${AUG_IMG:-}" ]; then
  echo "Probing Augustus image with $CONTAINER_RUNTIME: $AUG_IMG"
  # Try to run a lightweight command to ensure image is usable
  if ! $CONTAINER_RUNTIME exec "$AUG_IMG" augustus --version >/dev/null 2>&1; then
    echo "WARN: Could not run 'augustus --version' in $AUG_IMG; continuing but container may be unusable." >&2
  fi
  # Detect config directory regardless of docker:// or local .sif
  CFG_PATH=$(detect_cfg "$AUG_IMG") || true
  if [ -n "$CFG_PATH" ]; then
    echo "Using AUGUSTUS_CONFIG_PATH=$CFG_PATH"
    export SINGULARITYENV_AUGUSTUS_CONFIG_PATH="$CFG_PATH"
    export APPTAINERENV_AUGUSTUS_CONFIG_PATH="$CFG_PATH"
  else
    echo "WARN: Could not detect AUGUSTUS config dir in container; will rely on defaults"
  fi
fi

# Run multi-chromosome test wrapper (chr2 + chr18)
snakemake -s run_multi_test.smk \
  --configfile config/test_overlay_augustus.yaml \
  --profile config/slurm/profile \
  --use-conda \
  --use-singularity \
  --singularity-args "-B $PWD" \
  --conda-frontend "$FRONTEND" \
  --jobs 100 \
  --printshellcmds \
  target_multi
