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

# Detect container runtime (expect modules already loaded in the shell)
if command -v apptainer >/dev/null 2>&1; then
  RUNTIME=apptainer
elif command -v singularity >/dev/null 2>&1; then
  RUNTIME=singularity
else
  echo "ERROR: apptainer/singularity not found in PATH." >&2
  echo "Please load one in this shell, e.g.: 'module load apptainer' and rerun." >&2
  exit 1
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

# Verify the binary works inside the SIF (may require network for first pull)
echo "Verifying Augustus in $SIF_PATH"
$RUNTIME exec "$SIF_PATH" augustus --version

# Create/refresh a stable symlink for easier config
ln -sfn "$(basename "$SIF_PATH")" resources/containers/augustus.sif
STABLE_SIF="resources/containers/augustus.sif"

# Rewrite config to point containers.augustus at the stable SIF path
awk -v path="$STABLE_SIF" '
  BEGIN{f=0}
  /^containers:/{f=1}
  f && /^[[:space:]]+augustus:/{ sub(/:.*/, ": \"" path "\""); f=0 }
  {print}
' "$CONFIG_FILE" > "$CONFIG_FILE.tmp"
mv "$CONFIG_FILE.tmp" "$CONFIG_FILE"

cat <<MSG
SUCCESS: Augustus container is ready: $SIF_PATH

Next steps:
- The overlay has been updated to use: $STABLE_SIF
- Submit your Slurm job (no network required on compute nodes).
MSG
