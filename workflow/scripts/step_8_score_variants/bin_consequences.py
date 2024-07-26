#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 15-12-2023
:Usage: see basic_annotation.py --help

Generate basic annotation columns.
These can be derived from just the vcf and the reference genome.

These functions were copied (with minor modifications) from Gross et al.:
From pCADD VEP processing script:
- count_gc_cpg
From mCADD "Mouse-feature_annotate.py"
- parse_dna_shape
- get_shape_score (previously part of a larger function)
"""

# Import dependencies
import gzip
import math
import sys
from argparse import ArgumentParser

import pandas
import pysam
import numpy as np

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="Path to input .tsv file containing all annotations",
                    type=str, required=True)
parser.add_argument("-o", "--output",
                    help="Output file (default consequence_bins.tsv)",
                    type=str, default="basic_annotation.tsv")
parser.add_argument("-a", "--annotation",
                    help="CADD score annotation file, bed.gz format",
                    type=str, required=True)
args = parser.parse_args()

if __name__ == '__main__':
    outfile = gzip.open(args.output, "wt") if args.output.endswith(".gz") \
        else open(args.output, "w")

    df = pandas.read_csv(args.input, sep='\t', na_values=['-'])

    cadd_tabix = pysam.Tabixfile(args.annotation, 'r')

    ABBREVIATIONS = ["SG", "CS", "NS", "SN", "SL", "S", "U5", "U3", "R",
                     "IG", "NC", "I", "UP", "DN", "O"]

    bins = []
    for i in range(52):
        bins.append([0] * len(ABBREVIATIONS))
    # Annotate each variant
    for index, row in df.iterrows():
        pos = row["Pos"]
        found = False
        #print(f"Search: {pos}, {row['Ref']} {row['Alt']}")
        for elem in cadd_tabix.fetch(row["#Chrom"], pos, pos + 1):
            #print(elem)
            elem = elem.rstrip("\n").split('\t')
            if int(elem[1]) == pos and row["Ref"] == elem[2] and \
                    row["Alt"] == elem[3]:
                score = float(elem[5])
                consequence = row["Consequence"]
                #print(f"{pos}: {score}")
                bin = math.floor(score)
                if bin > 50:
                    bin = 51
                bins[bin][ABBREVIATIONS.index(consequence)] += 1
                found = True
                break
        if not found:
            sys.exit(f"No score found for: {pos}, {row['Ref']} {row['Alt']}")
    for bin, counts in enumerate(bins):
        outfile.write(str(bin) + "\t" +
                      ",".join([str(count) for count in counts]) + "\n")

    outfile.close()