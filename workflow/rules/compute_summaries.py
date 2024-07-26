#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 05-10-18

This script takes the fully computed and location sorted PHRED scores file
and computes per site summary statistics.

:Edited by: Job van Schipstal
:Date: 25-10-2023
:Usage: see compute_CADD_score_summary.py --help
"""

import sys
from argparse import ArgumentParser
import gzip

import numpy as np

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--infile", type=str, required=True,
                    help="sorted PHRED-score file. First line is ignored.")
parser.add_argument("-o", "--outmask", type=str,
                    help="PHRED-score containing file pattern, should "
                         "contain a mask TYPE, which will be replaced with "
                         "Median, Max, Min and Std for the different "
                         "operations that will be performed. "
                         "(Default: TYPE_PHRED_score)",
                    default="TYPE_PHRED_score.tsv")
args = parser.parse_args()

# Check that infile has TYPE mask
if "TYPE" not in args.infile:
    sys.exit(f"Mask sequence 'TYPE' not found in infile {args.infile}")

I_THRESHOLD = 3  # the 3 potential variants per position

# Open infile and skip header line
if args.infile.endswith("gz"):
    infile = gzip.open(args.infile, "rt")
else:
    infile = open(args.infile, "r")
infile_enumerator = enumerate(infile)
infile_enumerator.next()

# Open an outfile for each operation
names = ["Median", "Max", "Min", "Std"]
ops = [np.median, np.max, np.min, np.std]
operations = []
for name, operation in zip(names, ops):
    outfile = open(args.infile.replace("TYPE", name), "w")
    outfile.write(f"#Chrom\tPos\t{name}-PHRED\n")
    operations.append([name, operation, outfile])

# Loop through file and perform summary operations
Phred_store = []
for i, lines in infile_enumerator:
    line = lines.strip().split("\t")
    if (i % I_THRESHOLD == 0) and (i != 0):
        Phred_store.append(float(line[-1]))
        pos = line[0:2]
        pos = "\t".join(pos)
        [file.write(pos + "\t%s" % operation(Phred_store) + "\n")
         for name, operation, file in operations]
        Phred_store = []

    else:
        Phred_store.append(float(line[-1]))

# Close files
infile.close()
[file.close() for name, operation, file in operations]

"""
Test Input/Output:
X       1       C       A       0.5366852122925777      2.22782
X       1       C       G       0.5619843713968127      3.03890
X       1       C       T       0.5669613494334884      3.22649
X       2       A       C       0.5623888061466943      3.05383
X       2       A       G       0.5688997542352155      3.30157
X       2       A       T       0.5491479394495927      2.59782
X       3       G       A       0.48814700982315257     1.18826
X       3       G       C       0.6202772601258822      5.72309
X       3       G       T       0.5473617064825842      2.54151


#Chrom  Pos     Max-PHRED
X       1       3.22649
X       2       3.30157
X       3       5.72309

#Chrom  Pos     Min-PHRED
X       1       2.22782
X       2       2.59782
X       3       1.18826

#Chrom  Pos     Median-PHRED
X       1       3.0389
X       2       3.05383
X       3       2.54151
"""
