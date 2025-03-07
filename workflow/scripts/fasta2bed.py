#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

NUCLEOTIDES = {"A", "T", "C", "G"}

def process_fasta(infile: str) -> None:
    chrom = ""
    pos = -1
    start = -1
    in_sequence_region = False

    with open(infile, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if in_sequence_region:  # last sequence region from previous chrom
                    print(f"{chrom}\t{start}\t{pos}")
                    start = -1  # not needed actually
                    in_sequence_region = False
                pos = 0
                chrom = line.split(" ")[0].replace(">", "").strip()
                chrom = chrom.split(".")[-1]
            else:
                for c in line.strip():
                    if not in_sequence_region and c in NUCLEOTIDES:
                        in_sequence_region = True
                        start = pos
                    elif in_sequence_region and c not in NUCLEOTIDES:
                        in_sequence_region = False
                        print(f"{chrom}\t{start}\t{pos}")

                    pos += 1

    if in_sequence_region:  # last sequence region in last chrom
        print(f'{chrom}\t{start}\t{pos}')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: fasta2bed.py <input_fasta_file>")
        sys.exit(1)

    infile = sys.argv[1]
    try:
        process_fasta(infile)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)