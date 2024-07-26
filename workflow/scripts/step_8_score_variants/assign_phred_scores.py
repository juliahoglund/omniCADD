#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 21-09-18
This script takes the sorted (in descending order) RAW-scores and assigns to
them a PHRED-like score.

Percentage of variants within a certain score range:
0-10 : 90%
0-20 : 99%
0-30 : 99.9%
0-40 : 99.99%
0-50 : 99.999%
0-60 : 99.9999%
0-70 : 99.99999%
0-80 : 99.999999%
0-90 : 99.9999999%
0-100 : 100%
-10*math.log10(i/total)

:Edited by: Job van Schipstal
:Date: 25-10-2023
:Usage: see assign_phred_score.py --help

Modified to skip #header lines
and to write the scores to separate files per chromosome, easing filtering.
"""

import sys
import gzip
import math
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--infile", dest="raw",
                    help="input csv file, Chrom, Pos, Ref, Alt, RAW Score",
                    type=str, required=True)
parser.add_argument("-o", "--outmask",
                    help="output tsv file mask, with PHRED-like score added. "
                         "The mask CHROM is replaced with the actual "
                         "chromosome",
                    type=str, required=True)
parser.add_argument("--chroms",
                    help="list of chromosome names, multiple arguments "
                         "allowed. All chromosomes in the input file should "
                         "be passed",
                    type=str, required=True, nargs="+")
parser.add_argument("-c", "--counts", dest="line_counts",
                    help="int, total variant count", type=int, required=True)
args = parser.parse_args()


def open_file(mask, chrom):
    name = mask.replace("CHROM", chrom)
    if name.endswith(".gz"):
        phred_out = gzip.open(name, "wt")
    else:
        phred_out = open(name, "w")
    phred_out.write("#Chrom\tPos\tRef\tAlt\tRAW\tPHRED\n")
    return phred_out


# To ensure a float from the calculation, the PHRED score is not a whole number
line_counts = float(args.line_counts)

if args.raw.endswith(".gz"):
    raw_in = gzip.open(args.raw, "rt")
else:
    raw_in = open(args.raw, "r")

# Open separate output files for each chromosome
out_chrom = {}
for chromosome in args.chroms:
    out_chrom[chromosome] = open_file(args.outmask, chromosome)

# Loop though the input file, score all variants
idx = 0
for line in raw_in:
    if line.startswith("#"):
        continue
    idx += 1
    line = line.strip().split(",")
    score = -10 * math.log(idx / line_counts, 10)
    out_chrom[line[0]].write("\t".join(line) + f"\t{score:.8f}\n")

# Close files, good practice
raw_in.close()
[file.close() for file in out_chrom.values()]

# Verify that the amount of variants is as expected
if idx != line_counts:
    sys.exit(f"Ended at index: {idx} but expected {line_counts:.0f} variants!")
