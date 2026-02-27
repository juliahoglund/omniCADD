#!/usr/bin/env python3
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
from typing import Dict, TextIO

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--infile", dest="raw",
                    help="input csv file, Chrom, Pos, Ref, Alt, RAW Score",
# ...existing code...
