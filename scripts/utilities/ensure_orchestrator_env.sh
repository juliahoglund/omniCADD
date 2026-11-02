#!/bin/bash
# Create (if missing) and activate the Snakemake orchestrator venv, repairing it if
# snakemake or the pinned SLURM plugin are missing. Cluster independent: load whatever
# gives you python3 >= 3.11 first (module load python, or the system python).
#
# The interpreter is taken from $PYTHON (default python3) and must be >= 3.11 (Snakemake 9 needs it), e.g.
#   PYTHON=python3.13 source scripts/utilities/ensure_orchestrator_env.sh
# Usage (from a start script):  source scripts/utilities/ensure_orchestrator_env.sh [venv_dir]
# The venv defaults to <repo>/.venv; requirements live in workflow/envs/orchestrator-requirements.txt.
_repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
_venv="${1:-$_repo/.venv}"
_req="$_repo/workflow/envs/orchestrator-requirements.txt"

if [ ! -x "$_venv/bin/python" ]; then
    _py="${PYTHON:-python3}"
    if ! "$_py" -c 'import sys; sys.exit(sys.version_info < (3, 11))'; then
        echo "ERROR: $_py is $("$_py" --version 2>&1); Snakemake 9 needs Python >= 3.11. Set PYTHON=python3.x or load a newer python." >&2
        return 1
    fi
    echo "Creating orchestrator venv in $_venv with $_py"
    "$_py" -m venv "$_venv" || return 1
fi
source "$_venv/bin/activate"

_want_plugin=$(grep -E '^snakemake-executor-plugin-slurm==' "$_req" | cut -d= -f3)
# pip show exits non-zero when the package is missing; under the caller's `set -e -o pipefail` that would
# abort the whole start script, so tolerate it.
_have_plugin=$(pip show snakemake-executor-plugin-slurm 2>/dev/null | awk '/^Version:/{print $2}' || true)
if ! python -c "import snakemake" 2>/dev/null || [ "$_have_plugin" != "$_want_plugin" ]; then
    echo "Installing/repairing orchestrator packages from $_req"
    pip install --upgrade pip && pip install -r "$_req" || return 1
fi
echo "Using Snakemake from: $(command -v snakemake)"
