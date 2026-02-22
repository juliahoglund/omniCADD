#!/bin/bash
#SBATCH -A naiss2025-5-78
#SBATCH -J omniCADD_production
#SBATCH -t 48:00:00
#SBATCH -p main
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

# Production run: Full chromosome set with Augustus predictions

set -euo pipefail

# Load modules
ml PDC bioinfo-tools
ml snakemake/7.25.0 || true
ml apptainer/latest || true

# Select runtime binary
if command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
elif command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
else
    echo "ERROR: No container runtime found (singularity/apptainer)" >&2
    exit 1
fi
export APPTAINER_BIND="/cfs,/tmp"
export SINGULARITY_BIND="/cfs,/tmp"

# Check Augustus container
echo "Checking Augustus container: docker://juliahoglund/augustus:bioconda-3x"
$CONTAINER_CMD exec docker://juliahoglund/augustus:bioconda-3x augustus --version 2>&1 | grep -i "AUGUSTUS" && echo "[AUG] version rc= $?" || echo "[AUG] WARN: version check failed"
$CONTAINER_CMD exec docker://juliahoglund/augustus:bioconda-3x bash -c 'ls -d ${AUGUSTUS_CONFIG_PATH:-/opt/conda/envs/augustus/share/augustus/config} 2>/dev/null' && echo "[AUG] config: $($CONTAINER_CMD exec docker://juliahoglund/augustus:bioconda-3x bash -c 'echo ${AUGUSTUS_CONFIG_PATH:-/opt/conda/envs/augustus/share/augustus/config}')" || echo "[AUG] WARN: config not found"

# Check SNPEff container
echo "Checking SNPEff container: docker://juliahoglund/augustus:bioconda-3x"
$CONTAINER_CMD exec docker://juliahoglund/augustus:bioconda-3x snpEff -version 2>&1 | head -1 && echo "[SNP] version rc= $?" || echo "[SNP] WARN: version check failed"
$CONTAINER_CMD exec docker://juliahoglund/augustus:bioconda-3x java -version 2>&1 | head -1 && echo "[SNP] java rc= $?" || echo "[SNP] WARN: java not found"

if [ -f ".tests/snpEff/snpEff.config" ]; then
    echo "[SNP] config exists: True -> .tests/snpEff/snpEff.config"
else
    echo "[SNP] WARN: config missing: .tests/snpEff/snpEff.config"
fi

# Probe Augustus image with singularity
echo "Probing Augustus image with singularity: docker://juliahoglund/augustus:bioconda-3x"
$CONTAINER_CMD exec docker://juliahoglund/augustus:bioconda-3x bash -c '
  echo "Using AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH:-/opt/conda/envs/augustus/share/augustus/config}"
'

# Run snakemake with production config
snakemake -s run_production.smk \
  --configfile config/config_enhanced.yaml \
  --config-file config/production_augustus_full.yaml \
  --cores 16 \
  --use-singularity \
  --singularity-args "--bind /cfs,/tmp" \
  --use-conda \
  --conda-prefix "$PWD/.snakemake/conda" \
  --rerun-incomplete \
  --printshellcmds \
  --keep-going \
  --scheduler greedy
