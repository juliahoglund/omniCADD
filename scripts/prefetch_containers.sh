#!/usr/bin/env bash
set -euo pipefail

# Prefetch container images referenced in a config overlay.
# Usage:
#   ./scripts/prefetch_containers.sh [CONFIG_YAML]
# Default CONFIG_YAML: config/test_overlay_augustus.yaml
#
# This script pulls the Augustus container to a local SIF so compute nodes
# without outbound network can run with --use-singularity.

CONFIG_FILE=${1:-config/test_overlay_augustus.yaml}
if [ ! -f "$CONFIG_FILE" ]; then
  echo "Config file not found: $CONFIG_FILE" >&2
  exit 1
fi

# Load container runtime
module purge >/dev/null 2>&1 || true
module load apptainer >/dev/null 2>&1 || module load singularity >/dev/null 2>&1 || true
if ! command -v apptainer >/dev/null 2>&1 && ! command -v singularity >/dev/null 2>&1; then
  echo "ERROR: apptainer/singularity not available. Load the appropriate module (e.g., 'module load apptainer')." >&2
  exit 1
fi

if command -v apptainer >/dev/null 2>&1; then
  RUNTIME=apptainer
else
  RUNTIME=singularity
fi

AUG_IMG=$(awk '/^containers:/{f=1} f && /augustus:/ {print; exit}' "$CONFIG_FILE" | sed -E 's/.*"(.*)".*/\1/')
if [ -z "${AUG_IMG:-}" ]; then
  echo "No containers.augustus found in $CONFIG_FILE" >&2
  exit 1
fi

mkdir -p resources/containers

# Derive a SIF target path based on the tag if coming from docker://
if [[ "$AUG_IMG" == docker://* ]]; then
  TAG=$(echo "$AUG_IMG" | sed -E 's#.*/augustus:##')
  SAFE_TAG=$(echo "$TAG" | tr '/:' '__')
  SIF_PATH="resources/containers/augustus_${SAFE_TAG}.sif"
  echo "Pulling $AUG_IMG -> $SIF_PATH"
  $RUNTIME pull "$SIF_PATH" "$AUG_IMG"
else
  # Already a local path
  SIF_PATH="$AUG_IMG"
fi

# Verify the binary works inside the SIF
echo "Verifying Augustus in $SIF_PATH"
$RUNTIME exec "$SIF_PATH" augustus --version

cat <<MSG
SUCCESS: Augustus container is ready: $SIF_PATH

Next steps:
- Update containers.augustus in $CONFIG_FILE to the local SIF path:

  containers:
    augustus: "$SIF_PATH"

- Then submit your Slurm job (no network required on compute nodes).
MSG
