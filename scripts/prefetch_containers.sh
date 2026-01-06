#!/usr/bin/env bash
set -euo pipefail

# Prefetch container images referenced in a config overlay.
# Usage:
#   ./scripts/prefetch_containers.sh [CONFIG_YAML]
# Default CONFIG_YAML: config/test_overlay_augustus.yaml
#
# This script pulls the Augustus (and optionally SNPEff) containers to local SIFs
# so compute nodes without outbound network can run with --use-singularity.

CONFIG_FILE=${1:-config/test_overlay_augustus.yaml}
if [ ! -f "$CONFIG_FILE" ]; then
  echo "Config file not found: $CONFIG_FILE" >&2
  exit 1
fi

# Detect container runtime (best-effort auto-load PDC + apptainer)
if command -v apptainer >/dev/null 2>&1; then
  RUNTIME=apptainer
elif command -v singularity >/dev/null 2>&1; then
  RUNTIME=singularity
else
  # Try to load via PDC modules non-interactively
  if command -v ml >/dev/null 2>&1; then
    ml PDC >/dev/null 2>&1 || true
    ml apptainer >/dev/null 2>&1 || true
  elif command -v module >/dev/null 2>&1; then
    module load PDC >/dev/null 2>&1 || true
    module load apptainer >/dev/null 2>&1 || true
  fi
  if command -v apptainer >/dev/null 2>&1; then
    RUNTIME=apptainer
  elif command -v singularity >/dev/null 2>&1; then
    RUNTIME=singularity
  else
    echo "ERROR: apptainer/singularity not found in PATH." >&2
    echo "Please load modules, e.g.: 'ml PDC; ml apptainer' then rerun." >&2
    exit 1
  fi
fi

extract_img() {
  local name="$1"
  awk -v n="$name" '
    /^containers:/ {f=1; next}
    f && $0 ~ "^[[:space:]]*" n ":" {print; exit}
  ' "$CONFIG_FILE" | sed -E 's/.*"(.*)".*/\1/'
}

AUG_IMG=$(extract_img augustus)
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
  if [ -f "$SIF_PATH" ]; then
    echo "Found existing SIF: $SIF_PATH (reusing)"
  else
    echo "Pulling $AUG_IMG -> $SIF_PATH"
    if ! $RUNTIME pull "$SIF_PATH" "$AUG_IMG"; then
    echo "WARN: Failed to pull configured tag: $AUG_IMG" >&2
    echo "Attempting to auto-discover a valid Quay tag..." >&2
    # Query Quay for tags and try a handful of likely candidates.
    CANDIDATES=$(curl -fsSL 'https://quay.io/api/v1/repository/biocontainers/augustus/tag/?limit=200' \
      | python3 - <<'PY'
import sys, json
data=json.load(sys.stdin)
tags=[t['name'] for t in data.get('tags', [])]
# Prefer 3.x, then any
pref=[t for t in tags if t.startswith('3.')]
pref+= [t for t in tags if not t.startswith('3.')]
# Print top 20 candidates
print('\n'.join(pref[:20]))
PY
    ) || true
    if [ -z "$CANDIDATES" ]; then
      echo "ERROR: Could not retrieve tags from Quay. Check network or try later." >&2
      exit 1
    fi
    FOUND=""
    for t in $CANDIDATES; do
      CAND_URI="docker://quay.io/biocontainers/augustus:$t"
      SAFE=$(echo "$t" | tr '/:' '__')
      CAND_SIF="resources/containers/augustus_${SAFE}.sif"
      echo "Trying $CAND_URI"
      if $RUNTIME pull "$CAND_SIF" "$CAND_URI"; then
        echo "SUCCESS: Pulled $CAND_URI"
        SIF_PATH="$CAND_SIF"
        FOUND="$CAND_URI"
        break
      fi
    done
    if [ -z "$FOUND" ]; then
      echo "ERROR: Failed to pull any tested Augustus tags from Quay." >&2
      exit 1
    fi
    fi
  fi
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

# Optionally prefetch SNPEff if present
SNPEFF_IMG=$(extract_img snpeff || true)
if [ -n "${SNPEFF_IMG:-}" ]; then
  echo "\nDetected SNPEff container: $SNPEFF_IMG"
  # If SNPEff image equals Augustus image, reuse the same SIF and symlink
  if [ "$SNPEFF_IMG" = "$AUG_IMG" ] || [ "$SNPEFF_IMG" = "$STABLE_SIF" ]; then
    echo "SNPEff shares Augustus image; reusing $STABLE_SIF"
    ln -sfn "$(basename "$SIF_PATH")" resources/containers/snpeff.sif
    awk -v path="resources/containers/augustus.sif" '
      BEGIN{f=0}
      /^containers:/{f=1}
      f && /^[[:space:]]+snpeff:/{ sub(/:.*/, ": \"" path "\""); f=0 }
      {print}
    ' "$CONFIG_FILE" > "$CONFIG_FILE.tmp" && mv "$CONFIG_FILE.tmp" "$CONFIG_FILE"
  else
    SIF_SNP="resources/containers/snpeff.sif"
    if [[ "$SNPEFF_IMG" == docker://* ]]; then
      echo "Pulling $SNPEFF_IMG -> $SIF_SNP"
      $RUNTIME pull --force "$SIF_SNP" "$SNPEFF_IMG" || echo "WARN: Failed to pull SNPEff image; will rely on login-node cache at runtime" >&2
    else
      SIF_SNP="$SNPEFF_IMG"
    fi
    if [ -f "$SIF_SNP" ]; then
      echo "SNPEff SIF ready: $SIF_SNP"
      # Rewrite overlay to use stable SIF path
      awk -v path="resources/containers/snpeff.sif" '
        BEGIN{f=0}
        /^containers:/{f=1}
        f && /^[[:space:]]+snpeff:/{ sub(/:.*/, ": \"" path "\""); f=0 }
        {print}
      ' "$CONFIG_FILE" > "$CONFIG_FILE.tmp" && mv "$CONFIG_FILE.tmp" "$CONFIG_FILE"
      ln -sfn "$(basename "$SIF_SNP")" resources/containers/snpeff.sif
    fi
  fi
fi
