#!/bin/bash
# Aggregate Snakemake benchmark files into a summary
# Usage: ./aggregate_benchmarks.sh [benchmark_directory]

BENCH_DIR="${1:-logs/benchmarks}"

if [ ! -d "$BENCH_DIR" ]; then
    echo "Error: Benchmark directory '$BENCH_DIR' not found"
    echo "Note: Not all rules may have benchmarking enabled"
    exit 1
fi

echo "==================== BENCHMARK SUMMARY ===================="
echo ""
echo "Aggregating benchmark files from: $BENCH_DIR"
echo ""


echo "=== Resource Usage by Rule ==="
printf "%-30s %6s %12s %12s %10s %10s\n" "Rule" "Jobs" "Avg Runtime" "Total Runtime" "Avg RAM" "Max RAM"
echo "------------------------------------------------------------------------------------"

find "$BENCH_DIR" -name "*.tsv" -type f 2>/dev/null | \
while read -r bench_file; do
    # Extract rule name from filename (e.g., run_genome_vep_chr18_9.tsv -> run_genome_vep)
    basename "$bench_file" .tsv
done | \
sed 's/_chr[0-9X]*_.*//' | \
sort | uniq -c | \
while read -r count rule; do
    # Aggregate stats for this rule
    total_time=0
    total_mem=0
    max_mem=0
    job_count=0
    
    find "$BENCH_DIR" -name "${rule}*.tsv" -type f 2>/dev/null | \
    while read -r bench_file; do
        # first line is header
        tail -1 "$bench_file" 2>/dev/null
    done | \
    awk -v rule="$rule" '{
        # Column 1 is seconds
        total_time += $1
        # Column 3 is max_rss (RAM in MB)
        total_mem += $3
        if ($3 > max_mem) max_mem = $3
        job_count++
    }
    END {
        if (job_count > 0) {
            avg_time = total_time / job_count
            avg_mem = total_mem / job_count
            printf "%-30s %6d %12s %12s %9.0fM %9.0fM\n", 
                rule, job_count,
                format_time(avg_time),
                format_time(total_time),
                avg_mem, max_mem
        }
    }
    
    function format_time(sec) {
        if (sec < 60) return sprintf("%ds", sec)
        if (sec < 3600) return sprintf("%dm %ds", int(sec/60), sec%60)
        return sprintf("%dh %dm", int(sec/3600), int((sec%3600)/60))
    }'
done

echo ""
echo "=== Top 10 Longest Individual Jobs ==="
printf "%-50s %12s %10s\n" "Job" "Runtime" "Peak RAM"
echo "------------------------------------------------------------------------------------"

find "$BENCH_DIR" -name "*.tsv" -type f 2>/dev/null | \
while read -r bench_file; do
    job_name=$(basename "$bench_file" .tsv)
    tail -1 "$bench_file" 2>/dev/null | awk -v job="$job_name" '{
        printf "%s\t%s\t%s\n", job, $1, $3
    }'
done | \
sort -t$'\t' -k2 -rn | head -10 | \
awk -F'\t' '{
    job = $1
    time = $2
    mem = $3
    
    # Format time
    if (time < 60) time_str = sprintf("%ds", time)
    else if (time < 3600) time_str = sprintf("%dm %ds", int(time/60), time%60)
    else time_str = sprintf("%dh %dm", int(time/3600), int((time%3600)/60))
    
    printf "%-50s %12s %9.0fM\n", substr(job, 1, 50), time_str, mem
}'

echo ""
echo "================================================================"
