#!/usr/bin/env python3
# -*- coding: ASCII -*-

"""
:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 15-08-18

This script takes fully annotated variant files for simulated and derived
variants respectively. The script encodes and imputes data of features.

:Edited by: Job van Schipstal, Julia Höglund
:Date: 16-10-2023; 2025-03-05
:Usage: see data_preparation.py --help

Modifications:
- Script is configured via tsv files, instead of hard-coding
  which columns should be processed in what way.
- Updated from optparse to argparse
- Modified to work on parts of the variants (e.g. only a single chromosome)
  of either derived,simulated or whole genome variants.
- Deriving mean functionality is implemented in a separate script,
  derive_means.py so that after that step completes all variants can be
  processed in parallel using this script.
- The interaction term generation is done with sparse floats
  instead of dense values, to dramatically reduce memory usage.
  (density was only 5% in testing). Afterwards the dataset can also be saved
  as an scipy sparse matrix, which is directly supported by scikit-learn.
"""

# Import dependencies
import sys
from argparse import ArgumentParser

import pandas
import ast
import numpy as np
from scipy.sparse import csr_matrix, save_npz

# ...existing code...
