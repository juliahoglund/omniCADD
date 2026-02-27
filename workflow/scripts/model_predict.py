#!/usr/bin/env python3
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 25-10-2023
:Usage: see test_model.py --help

Loads model chunks from npz files, merges chunks and then scales by standard
deviation but not mean centered, as this would break sparsity. Writes the
predicted probability for a variant to be of class 1, (proxy) deleterious.
Additionally, Chrom, Pos, Ref- and Alternative nucleotide are in the output.
"""

# Import dependencies
import sys
import pickle
from itertools import product
from argparse import ArgumentParser

import pandas
import numpy as np
from scipy.sparse import load_npz, vstack
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="test set variants, npz file(s), at least one should "
                         "be provided. It is expected that there are "
                         ".meta.csv.gz and columns.csv files present as well",
                    type=str, required=True, nargs="+")
parser.add_argument("--model",
                    help="Pickle file containing the scikit-learn model",
                    type=str, required=True)
parser.add_argument("--scaler",
                    help="Pickle file containing the scikit-learn scaler",
                    type=str, required=True)
parser.add_argument("-o", "--output",
# ...existing code...
