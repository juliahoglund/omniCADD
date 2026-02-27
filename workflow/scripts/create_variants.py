#!/usr/bin/env python3
# -*- coding: ASCII -*-

"""
:Author: Christian Geoss
:Contact: c.gross@tudelft.nl
:Date: 06.05.2019

:Edited by: Job van Schipstal
:Date: 24-10-2023
:Usage: see create_variants.py --help

This script iterates over the reference sequence and generates all possible
SNPs for the chromosome of interest. These are stored in .vcf files of a
user specified size at the desired location.
"""

import re
import sys
from argparse import ArgumentParser
from typing import Iterator

from Bio import SeqIO

parser = ArgumentParser(__doc__)
parser.add_argument("-o", "--out-path", type=str,
                    help="Path to output folder, "
                         "(default: the working directory)", default="")
parser.add_argument("-s", "--size", type=int,
                    help="Number of sites for which variants are created. "
                         "The total number of lines in the output is 3 times "
                         "the specified number, for the 3 alternative alleles "
                         "per site. (default: 1000000)", default=1000000)
parser.add_argument("-r", "--reference", type=str,
                    help="Path to reference, fasta format", required=True)
parser.add_argument("-p", "--start-pos", type=int,
                    help="The position at which to start generating variants, "
                         "0-based (default: 0)", default=0)
parser.add_argument("-c", "--chrom", type=str,
                    help="The chromosome for which to generate variants, "
# ...existing code...
