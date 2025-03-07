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
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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

# Constants
ABBREVIATIONS = ["SG", "CS", "NS", "SN", "SL", "S", "U5", "U3", "R",
                 "IG", "NC", "I", "UP", "DN", "O"]
NUM_BINS = 52

def initialize_bins(num_bins, abbreviations):
    return [[0] * len(abbreviations) for _ in range(num_bins)]

def process_variant(row, cadd_tabix, bins, abbreviations):
    pos = row["Pos"]
    found = False
    for elem in cadd_tabix.fetch(row["#Chrom"], pos, pos + 1):
        elem = elem.rstrip("\n").split('\t')
        if int(elem[1]) == pos and row["Ref"] == elem[2] and row["Alt"] == elem[3]:
            score = float(elem[5])
            consequence = row["Consequence"]
            bin = min(math.floor(score), 51)
            bins[bin][abbreviations.index(consequence)] += 1
            found = True
            break
    if not found:
        logging.error(f"No score found for: {pos}, {row['Ref']} {row['Alt']}")
        sys.exit(1)

def write_output(outfile, bins):
    for bin, counts in enumerate(bins):
        outfile.write(str(bin) + "\t" + ",".join([str(count) for count in counts]) + "\n")

if __name__ == '__main__':
    try:
        outfile = gzip.open(args.output, "wt") if args.output.endswith(".gz") else open(args.output, "w")
    except IOError as e:
        logging.error(f"Error opening output file: {e}")
        sys.exit(1)

    try:
        df = pandas.read_csv(args.input, sep='\t', na_values=['-'], low_memory=False)
    except Exception as e:
        logging.error(f"Error reading input file: {e}")
        sys.exit(1)

    try:
        cadd_tabix = pysam.Tabixfile(args.annotation, 'r')
    except Exception as e:
        logging.error(f"Error opening annotation file: {e}")
        sys.exit(1)

    bins = initialize_bins(NUM_BINS, ABBREVIATIONS)

    for index, row in df.iterrows():
        process_variant(row, cadd_tabix, bins, ABBREVIATIONS)

    write_output(outfile, bins)
    outfile.close()
    logging.info("Processing completed successfully.")