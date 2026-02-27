#!/usr/bin/env python3
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 7-12-2023
:Usage: see fold_data.py --help

Loads model chunks from npz files, and the y variable from the metadata csv.
Merges chunks and determines the folds.
Each fold is then written to a separate file.
"""

# Import dependencies
import sys
from argparse import ArgumentParser
from pathlib import Path

from joblib import Parallel, delayed
import sklearn
from sklearn.model_selection import StratifiedKFold
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="variants, npz file(s), at least one should be provided. "
                         "It is expected that there are .meta.csv.gz and columns.csv files present as well",
                    type=str,
                    required=True,
                    nargs="+")
parser.add_argument("-o", "--output",
                    help="output fold files, the amount of files provided equals the number of folds, "
                         "would not recommend going lower as 4",
                    type=str,
                    required=True,
                    nargs="+")
# ...existing code...
