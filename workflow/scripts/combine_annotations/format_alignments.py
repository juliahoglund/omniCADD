#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Linearize a per-chunk multi-block FASTA (from maf2fasta.pl) into one
continuous sequence per species, gap-filling species missing from
individual blocks, and write an index of genomic positions for the
species of interest.

:Modified by: Julia Höglund 2025-03-05
:Rewritten: 2026-08-03
"""

import sys


def parse_blocks(ffile):
    """
    Parse maf2fasta.pl output into a list of blocks.
    Each block is a dict: species_id -> (start, sequence).
    """
    blocks = []
    current = {}
    current_id = None
    current_start = None
    current_seq_parts = []

    def flush_record():
        if current_id is not None:
            current[current_id] = (current_start, "".join(current_seq_parts))

    with open(ffile, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if line == "=":
                flush_record()
                blocks.append(current)
                current = {}
                current_id = None
                current_start = None
                current_seq_parts = []
            elif line.startswith(">"):
                flush_record()
                fields = line[1:].split()
                current_id = fields[0].split(".")[0]
                current_start = int(fields[1])
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

    # in case the file doesn't end with a trailing "=" block separator
    if current:
        flush_record()
        blocks.append(current)

    return blocks


def main():
    if len(sys.argv) != 5:
        sys.exit("usage: format_alignments.py FASTA OUTPUT INDEXFILE SPECIES_OF_INTEREST")

    ffile, output, indexfile, species_of_interest = sys.argv[1:5]

    blocks = parse_blocks(ffile)

    # species order: species of interest first (matches prune_columns.py,
    # which prunes columns based on the FIRST record in the alignment),
    # then all others in first-appearance order
    species_order = [species_of_interest]
    seen = {species_of_interest}
    for block in blocks:
        for sp in block:
            if sp not in seen:
                seen.add(sp)
                species_order.append(sp)

    sequences = {sp: [] for sp in species_order}
    index_lines = []

    for block_index, block in enumerate(blocks):
        if not block:
            continue
        widths = {len(seq) for _, seq in block.values()}
        if len(widths) > 1:
            per_species = {sp: len(seq) for sp, (_, seq) in block.items()}
            sys.exit(
                f"ERROR: inconsistent sequence widths within block {block_index} "
                f"(of {len(blocks)} blocks): {per_species}"
            )
        width = widths.pop()

        for sp in species_order:
            if sp in block:
                sequences[sp].append(block[sp][1])
            else:
                sequences[sp].append("-" * width)

        if species_of_interest in block:
            start, _ = block[species_of_interest]
            index_lines.append(f"start: {start}, size: {width}\n")

    with open(output, "w") as fasta_out:
        for i, sp in enumerate(species_order):
            if i > 0:
                fasta_out.write("\n")
            fasta_out.write(">" + sp + "\n")
            fasta_out.write("".join(sequences[sp]))

    with open(indexfile, "w") as index_out:
        index_out.writelines(index_lines)


if __name__ == "__main__":
    main()
