#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

from data_helper import process_fasta

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: fasta2bed.py <input_fasta_file>")
        sys.exit(1)
    infile = sys.argv[1]
    regions = process_fasta(infile)
    for chrom, start, end in regions:
        print(f"{chrom}\t{start}\t{end}")