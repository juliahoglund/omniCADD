#!/bin/bash
# Test script to validate workflow with both standard and limited data configurations
# Usage: bash .tests/test_configs.sh

set -e  # Exit on error
set -o pipefail  # Pipelines fail if any command fails

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "================================================"
echo "  omniCADD Configuration Testing"
echo "================================================"
echo ""

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to run dry-run test
test_config() {
    local config_file=$1
    local config_name=$2
    
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo -e "${YELLOW}Testing: $config_name${NC}"
    echo "Config: $config_file"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    cd "$REPO_ROOT"
    
    # Copy test config to main config location
    cp "$config_file" config/config.yaml
    
    # Ensure other configs exist
    if [ ! -f config/config_refs.yaml ]; then
        cp config/config_refs.default.yaml config/config_refs.yaml
    fi
    if [ ! -f config/resources.yaml ]; then
        cp config/resources.default.yaml config/resources.yaml
    fi
    
    # Run dry-run
    echo ""
    echo "Running: snakemake --dry-run test_dag --configfile config/config.yaml"
    
    # Capture output and exit code
    set +e  # Temporarily disable exit on error to capture status
    snakemake --dry-run test_dag --configfile config/config.yaml > "/tmp/snakemake_${config_name}.log" 2>&1
    EXIT_CODE=$?
    set -e
    
    # Display output
    cat "/tmp/snakemake_${config_name}.log"
    
    # Check result
    if [ $EXIT_CODE -eq 0 ]; then
        echo ""
        echo -e "${GREEN}✓ $config_name: PASSED${NC}"
        return 0
    else
        echo ""
        echo -e "${RED}✗ $config_name: FAILED (exit code: $EXIT_CODE)${NC}"
        echo "See log: /tmp/snakemake_${config_name}.log"
        return 1
    fi
}

# Check if Snakemake is installed
if ! command -v snakemake &> /dev/null; then
    echo -e "${RED}Error: Snakemake not found${NC}"
    echo "Please install: conda install -c conda-forge -c bioconda snakemake"
    exit 1
fi

echo "Snakemake version: $(snakemake --version)"
echo ""

# Track results
PASSED=0
FAILED=0

# Test standard configuration
if test_config "$SCRIPT_DIR/config/config_standard.yaml" "standard"; then
    PASSED=$((PASSED + 1))
else
    FAILED=$((FAILED + 1))
fi

# Test limited configuration
if test_config "$SCRIPT_DIR/config/config_limited.yaml" "limited"; then
    PASSED=$((PASSED + 1))
else
    FAILED=$((FAILED + 1))
fi

# Summary
echo ""
echo "================================================"
echo "  Test Summary"
echo "================================================"
echo -e "Passed: ${GREEN}$PASSED${NC}"
echo -e "Failed: ${RED}$FAILED${NC}"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed${NC}"
    exit 1
fi
