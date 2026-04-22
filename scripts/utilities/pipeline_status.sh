#!/bin/bash
# Comprehensive pipeline status checker
# Usage: ./pipeline_status.sh [--timing] [--benchmarks]

SHOW_TIMING=0
SHOW_BENCHMARKS=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --timing|-t) SHOW_TIMING=1; shift ;;
        --benchmarks|-b) SHOW_BENCHMARKS=1; shift ;;
        --all|-a) SHOW_TIMING=1; SHOW_BENCHMARKS=1; shift ;;
        *) LOG_FILE="$1"; shift ;;
    esac
done

LOG_FILE="${LOG_FILE:-snakemake_run.log}"

echo "==================== PIPELINE STATUS ===================="
echo ""

# Basic progress from log
if [ -f "$LOG_FILE" ]; then
    echo "=== Job Statistics (from DAG) ==="
    grep -A 20 "^Job stats:" "$LOG_FILE" 2>/dev/null | grep -v "^--" | head -30
    
    echo ""
    echo "=== Progress Overview ==="
    COMPLETED=$(grep -c "Finished job" "$LOG_FILE" 2>/dev/null)
    TOTAL=$(grep "^total " "$LOG_FILE" | awk '{print $2}' | tail -1)
    
    if [ -n "$TOTAL" ] && [ "$TOTAL" -gt 0 ]; then
        PERCENT=$(echo "scale=1; ($COMPLETED * 100) / $TOTAL" | bc 2>/dev/null || echo "?")
        echo "Completed: $COMPLETED / $TOTAL jobs ($PERCENT%)"
    else
        echo "Completed: $COMPLETED jobs"
    fi
    
    echo ""
    echo "=== Recent Activity ==="
    echo "Last 5 completions:"
    grep "Finished job" "$LOG_FILE" | tail -5 | sed 's/.*Finished job /  /'
    
    ERRORS=$(grep -c "Error in rule" "$LOG_FILE" 2>/dev/null)
    if [ "$ERRORS" -gt 0 ]; then
        echo ""
        echo "⚠️  Errors detected: $ERRORS"
        echo "Recent errors:"
        grep "Error in rule" "$LOG_FILE" | tail -3 | sed 's/Error in rule /  /'
    fi
else
    echo "Log file not found: $LOG_FILE"
fi

# Show timing if requested
if [ "$SHOW_TIMING" -eq 1 ] && [ -f "$LOG_FILE" ]; then
    echo ""
    echo "========================================================="
    scripts/utilities/timing_summary.sh "$LOG_FILE"
fi

# Show benchmarks if requested  
if [ "$SHOW_BENCHMARKS" -eq 1 ] && [ -d "logs/benchmarks" ]; then
    echo ""
    echo "========================================================="
    scripts/utilities/aggregate_benchmarks.sh
fi

echo ""
echo "========================================================="
echo ""
echo "Tips:"
echo "  - Follow live: tail -f $LOG_FILE"
echo "  - Show timing: $0 --timing"
echo "  - Show benchmarks: $0 --benchmarks"
echo "  - Show all: $0 --all"
