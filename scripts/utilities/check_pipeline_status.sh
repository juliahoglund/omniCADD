#!/bin/bash
# Parse Snakemake log to show pipeline progress
# Usage: ./check_pipeline_status.sh [log_file]

LOG_FILE="${1:-snakemake_run.log}"

if [ ! -f "$LOG_FILE" ]; then
    echo "Error: Log file '$LOG_FILE' not found"
    echo "Usage: $0 [log_file]"
    exit 1
fi

echo "==================== PIPELINE STATUS ===================="
echo ""

# Extract job stats from the initial DAG planning
echo "=== Job Statistics (from DAG) ==="
grep -A 20 "^Job stats:" "$LOG_FILE" 2>/dev/null | grep -v "^--" | head -30

echo ""
echo "=== Progress Overview ==="
# Count completed jobs from log
COMPLETED=$(grep -c "Finished job" "$LOG_FILE" 2>/dev/null)
# Get total jobs from job stats if available
TOTAL=$(grep "^total " "$LOG_FILE" | awk '{print $2}' | tail -1)

if [ -n "$TOTAL" ] && [ "$TOTAL" -gt 0 ]; then
    PERCENT=$(echo "scale=1; ($COMPLETED * 100) / $TOTAL" | bc 2>/dev/null || echo "?")
    echo "Completed: $COMPLETED / $TOTAL jobs ($PERCENT%)"
else
    echo "Completed: $COMPLETED jobs"
fi

echo ""
echo "=== Currently Running/Submitted ==="
grep "rule " "$LOG_FILE" | grep -v "Finished job" | tail -10

echo ""
echo "=== Recent Completions ==="
grep "Finished job" "$LOG_FILE" | tail -10

echo ""
echo "=== Recent Errors/Failures ==="
grep -E "(Error|failed|FAILED)" "$LOG_FILE" | tail -5 || echo "No errors found"

echo ""
echo "=== Last Log Update ==="
ls -lh "$LOG_FILE" | awk '{print "Size: " $5 " | Modified: " $6 " " $7 " " $8}'

echo ""
echo "========================================================="
echo "Tip: Use 'tail -f $LOG_FILE' to follow progress in real-time"
