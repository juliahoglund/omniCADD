#!/bin/bash
# Parse Snakemake log for timing statistics
# Usage: ./timing_summary.sh [snakemake_run.log]

LOG_FILE="${1:-snakemake_run.log}"

if [ ! -f "$LOG_FILE" ]; then
    echo "Error: Log file '$LOG_FILE' not found"
    exit 1
fi

echo "==================== TIMING SUMMARY ===================="
echo ""

# Extract job completion times and calculate statistics per rule
echo "=== Average Runtime by Rule ==="
grep "Finished job" "$LOG_FILE" 2>/dev/null | \
    awk '{
        # Extract rule name and time
        for(i=1; i<=NF; i++) {
            if($i == "rule") rule=$(i+1)
        }
        # Time is usually at the end like "5h 30m 45s" or "45s" or "2m 30s"
        time_str = $(NF-1) " " $NF
        
        # Convert to seconds (simplified - handles s, m, h)
        seconds = 0
        if (match(time_str, /([0-9]+)h/, arr)) seconds += arr[1] * 3600
        if (match(time_str, /([0-9]+)m/, arr)) seconds += arr[1] * 60
        if (match(time_str, /([0-9]+)s/, arr)) seconds += arr[1]
        
        if (seconds > 0 && rule != "") {
            rule_times[rule] += seconds
            rule_counts[rule]++
            total_time += seconds
            total_jobs++
        }
    }
    END {
        # Sort by total time descending
        n = asorti(rule_times, sorted_rules, "@val_num_desc")
        
        printf "%-30s %8s %10s %12s %10s\n", "Rule", "Jobs", "Total", "Average", "% Time"
        printf "%s\n", "--------------------------------------------------------------------------------"
        
        for (i = 1; i <= n; i++) {
            rule = sorted_rules[i]
            total = rule_times[rule]
            count = rule_counts[rule]
            avg = total / count
            pct = (total / total_time) * 100
            
            printf "%-30s %8d %10s %12s %9.1f%%\n", 
                rule, count, 
                format_time(total), 
                format_time(avg), 
                pct
        }
        
        printf "%s\n", "--------------------------------------------------------------------------------"
        printf "%-30s %8d %10s\n", "TOTAL", total_jobs, format_time(total_time)
    }
    
    function format_time(sec) {
        if (sec < 60) return sprintf("%ds", sec)
        if (sec < 3600) return sprintf("%dm %ds", int(sec/60), sec%60)
        return sprintf("%dh %dm", int(sec/3600), int((sec%3600)/60))
    }'

echo ""
echo "=== Pipeline Duration ==="
START_TIME=$(grep -m1 "Building DAG" "$LOG_FILE" | awk '{print $1, $2}')
END_TIME=$(grep "of [0-9]* steps.*done" "$LOG_FILE" | tail -1 | awk '{print $1, $2}')

if [ -n "$START_TIME" ] && [ -n "$END_TIME" ]; then
    echo "Started:  $START_TIME"
    echo "Finished: $END_TIME"
else
    echo "Pipeline still running or log incomplete"
fi

echo ""
echo "========================================================="
