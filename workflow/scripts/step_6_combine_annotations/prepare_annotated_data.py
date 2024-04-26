#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 15-08-18

This script takes fully annotated variant files for simulated and derived
variants respectively. The script encodes and imputes data of features.

:Edited by: Job van Schipstal
:Date: 16-10-2023
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

"""
:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 15-08-18

This script takes fully annotated variant files for simulated and derived
variants respectively. The script encodes and imputes data of features.

Instead of reading in and manipulating both files at the same time it will
be performed one at the time, but the simulated data first because mean
imputation is used from the simulated data. It check automatically if
imputation file is present in current directory, if not and if derived is
specified then an error is thrown.

:Edited by: Seyan Hu
:Date: 13-12-2022
:Usage: python data_encoding_mv_handling.py -i <Merged annotation file>
(add '-d' if this is performed on derived variants)

"""

# Import dependencies
import sys
from argparse import ArgumentParser

import pandas
import ast
import numpy as np
from scipy.sparse import csr_matrix, save_npz

# Create and process cli
parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="annotated infile containing annotated variants to process", 
                    type=str, 
                    required=True)

parser.add_argument("-c", "--csv",
                    help="Output file, processed dataset in csv format (default None, not written)",
                    type=str, 
                    default="None")

parser.add_argument("-n", "--npz",
                    help="NPZ format output, .npz will be appended if not "
                         "already present. Metadata is stored seperately in "
                         ".meta.csv.gz and the columns in columns.csv (default"
                         " None, not written)",
                    type=str, 
                    default="None")

parser.add_argument("--processing-config",
                    help="Configuration tsv file indicating how to dataset should be processed",
                    type=str, 
                    required=True)

parser.add_argument("--interaction-config",
                    help="Configuration tsv file indicating which "
                         "interaction terms should be generated in the dataset",
                    type=str, 
                    required=True)

parser.add_argument("--imputation-dict",
                    help="Dictionary file to read/write imputation values to.",
                    type=str, 
                    default="impute_dict.txt")

parser.add_argument("-d", "--derived",
                    help="Is this derived data? Some columns are renamed"
                         "(default: False)",
                    default=False, 
                    action="store_true")

parser.add_argument("-y", "--y-value",
                    help="Y_value of data, added in new column y, must be float (default None, column not added)",
                    type=float, 
                    default=None)

args = parser.parse_args()

if args.csv == "None" and args.npz == "None":
    sys.stderr.write("WARNING: No output filetype is selected, performing "
                     "trial data preparation but not saving results.")


def load_tsv_configuration(file: str) -> dict:
    """
    Loads configuration from tsv table file.
    The first non # or empty line is read as column names

    :param file: str, filename to load configuration from
    :return: dict key is entry label and value is dict,
    with key column label and the value of that column
    """
    file_h = open(file, "r")
    elements = None
    samples = {}
    for line in file_h:
        line = line.strip()
        parts = line.split("\t")
        if line.startswith("#") or len(line) == 0:
            continue
        if not elements:
            elements = parts[1:]
            continue
        samples[parts[0]] = dict(
            [(key, value) for key, value in zip(elements, parts[1:])])
    return samples

def class_encoded_check(classlabel, selection, data):
    """
    Checks whether one-hot encoding resulted in all variants, since missing
    columns would pose an issue during when running the model.
    :param classlabel: str, label of original annotation, used as prefix
    :param selection: Iterable, columns to check existence of
    :param data: dataframe to check encoding in
    :return: dataframe that has been encoded
    """
    for clazz in selection:
        col = classlabel + '_' + clazz
        try:
            data[col]
        except KeyError:
            data[col] = np.zeros(data.shape[0], dtype=float)
    return data


 