#!/usr/bin/env python

import sys

chrom = ""
pos = -1
start = -1
in_sequence_region = False
infile = sys.argv[1]
with open(sys.argv[1], "r") as fh:
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
                if not in_sequence_region and c in ["A", "T", "C", "G"]:
                    in_sequence_region = True
                    start = pos
                elif in_sequence_region and c not in ["A", "T", "C", "G"]:
                    in_sequence_region = False
                    print(f"{chrom}\t{start}\t{pos}")

                pos += 1

if in_sequence_region:  # last sequence region in last chrom
    print(f"{chrom}\t{start}\t{pos}")