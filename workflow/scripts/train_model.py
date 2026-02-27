#!/usr/bin/env python3
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 11-12-2023
:Usage: see train_test_model.py --help

Loads model chunks from npz files, and the y variable from the metadata csv.
Merges chunks and splits them into a train and test set, the data is then
scaled by standard deviation but not mean centered, as this would break
sparsity. Trains a logistic regression model on the train set and generates
basic statistics for the performance on the test set.

Will save the train or test dataset, the model or the scaler,
or a roc-auc plot if these arguments are given.
"""

# Import dependencies
import sys
import os
import pickle
import logging
from argparse import ArgumentParser
from pathlib import Path

import pandas
from joblib import Parallel, delayed

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, confusion_matrix, \
    classification_report

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = ArgumentParser(description=__doc__)
parser.add_argument("--train",
# ...existing code...
