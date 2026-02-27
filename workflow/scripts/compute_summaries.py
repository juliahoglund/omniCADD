#!/usr/bin/env python3
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

# ...existing code...
