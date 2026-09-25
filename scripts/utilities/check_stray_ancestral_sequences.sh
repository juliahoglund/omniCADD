#!/bin/bash
# Post-hoc check for the mark_ancestor.py fix (stray "ancestral_sequences"
# pseudo-species leaking through for blocks where the target ancestor
# wasn't found/matched - see workflow/scripts/ancestral_generation/mark_ancestor.py).
# Confirms no "s ancestral_sequences.*" sequence row remains in the final
# sorted alignments after a full rerun.
#
# Usage: run from the pipeline working directory (e.g.
# /cfs/klemming/.../omniCADD/restructure) after the rerun completes:
#   bash scripts/utilities/check_stray_ancestral_sequences.sh

set -euo pipefail

fail=0
total_stray=0

for f in results/alignment/sorted/chr*.maf.lz4; do
    [ -e "$f" ] || continue
    count=$(lz4 -dc "$f" | grep -c "^s ancestral_sequences\." || true)
    total_stray=$((total_stray + count))
    if [ "$count" -gt 0 ]; then
        echo "FAIL: $f has $count stray ancestral_sequences record(s)"
        fail=1
    else
        echo "OK:   $f clean"
    fi
done

echo
if [ "$fail" -eq 0 ]; then
    echo "All sorted alignments clean - no stray ancestral_sequences records found."
else
    echo "Found $total_stray stray record(s) total across the files listed above."
    echo "This means mark_ancestor.py's fix isn't fully reflected here - check"
    echo "whether these chromosomes actually reran after the fix was pulled"
    echo "(rerun-triggers: input should have forced it), or whether the fix"
    echo "itself has a remaining edge case."
fi

exit "$fail"
