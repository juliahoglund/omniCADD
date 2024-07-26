#!/usr/bin/env python
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
                    help="File to write scored variants to, they will be in "
                         "uncompressed csv fomat (default: scored.csv)",
                    type=str, default="scored.csv")
parser.add_argument("--sort",
                    help="Sort scores in descending order",
                    action="store_true")
parser.add_argument("--no-header",
                    help="Do not write header in the output file",
                    action="store_true")
args = parser.parse_args()

NUCLEOTIDES = ["A", "C", "T", "G"]


def load_npz_with_meta(files) -> tuple:
    """
    Loads data in npz format to scipy sparse matrix and metadata in csv
    format to a pandas dataframe, column list is also returned and it is
    verified that this is the same for each file.
    :param files: Iterable of str, files to load the data from.
    :return: tuple, scipy sparse matrix, pandas dataframe, list of column names
    """
    data = []
    meta = []
    cols = None
    for file in files:
        sys.stderr.write(f"## Loading data from {file.split('/')[-1]}\n")
        # Load cols and verify they match for the different files
        new_cols = open(file + ".columns.csv").read().strip().split(",")
        if not cols:
            cols = new_cols
        elif not cols == new_cols:
            for i in range(min(len(cols), len(new_cols))):
                if cols[i] != new_cols[i]:
                    sys.stderr.write(f"{i}: {cols[i]} vs {new_cols[i]}\t")
            sys.exit(f"\nColumns for file {file} differ from those of the "
                     f"previously read infile(s), "
                     f"lengths: {len(new_cols)} vs {len(cols)}\n" +
                     " ".join(new_cols) + "\n" + " ".join(cols) + "\n")

        # Load metadata via pandas df and the actual data via scipy sparse.
        meta.append(pandas.read_csv(file + ".meta.csv.gz", sep=',',
                                    na_values=['-'], header=0))
        data.append(load_npz(file))
    sys.stderr.write("## Merging data parts\n")
    # CSR is row based, stacking it by row should be efficient, logistic
    # regression from scikit-learn also expects a csr format sparse matrix
    # so we want it to remain of that type.
    data = vstack(data, format="csr")
    return data, pandas.concat(meta, axis=0), cols


def get_pickle(file, clazz):
    """
    Loads an object from file and checks if it is an instance of clazz
    :param file: file to read object from
    :param clazz: clazz that is expected to be returned
    :return: instance of clazz, loaded from file
    """
    with open(file, "rb") as file_h:
        ins = pickle.load(file_h)
    if not isinstance(ins, clazz):
        sys.exit(f"Error loading instance from picle file, "
                 f"expected an instance of {clazz} but got {type(ins)}")
    return ins


def get_df(desired, all_cols, d_mat):
    """
    Returns desired columns from a scipy data matrix as a dense pandas df.

    :param desired: list of str, desired columns
    :param all_cols: list of str, column names of data_matrix
    :param d_mat: Scipy sparse matrix  to extract columns from
    :return: Pandas dataframe with desired columns
    """
    id_cols = {}
    for col in desired:
        if col not in all_cols:
            sys.exit(f"Missing column {col} in data_matrix, "
                     f"required for this analysis")
        id_cols[all_cols.index(col)] = col
    id_cols = dict(sorted(id_cols.items()))
    data_sel = d_mat[:, list(id_cols.keys())]
    print(list(id_cols.values()))
    return pandas.DataFrame(data_sel.todense(), columns=list(id_cols.values()))


def decode(n_tuple, prefix):
    """
    Retrieves nucleotide from one-hot-encoded nt columns.
    :param n_tuple: named tuple, row to find column in
    :param prefix: str, prefix of column name
    :return: str nucleotide or numpy's Nan if none found
    """
    for nt in NUCLEOTIDES:
        if getattr(n_tuple, prefix + nt) != 0:
            return nt
    sys.stderr.write(f"Unable to find nt in {tuple}\n")
    return np.NaN


# Load dataset
data_matrix, meta_data, columns = load_npz_with_meta(args.input)

# Extract needed metadata
chromosome = meta_data[["Chrom"]].to_numpy().flatten()
pos = meta_data[["Pos"]].to_numpy().flatten()
del meta_data

model = get_pickle(args.model, LogisticRegression)
scaler = get_pickle(args.scaler, StandardScaler)

# Scale and predict class
# sys.stderr.write("## Scaling dataset and predicting class\n")
scaled_data = scaler.transform(data_matrix)
del data_matrix
y_pred_prob = model.predict_proba(scaled_data)[::, 1]

# Get ref and alt cols since these are desired in the output
cols = [pref + nt for pref, nt in product(["Ref_", "Alt_"], NUCLEOTIDES)]
sub_df = get_df(cols, columns, scaled_data)
del scaled_data
ref = [decode(row, "Ref_") for row in sub_df.itertuples(index=False)]
alt = [decode(row, "Alt_") for row in sub_df.itertuples(index=False)]
del sub_df

df = pandas.DataFrame.from_dict({"#Chrom": chromosome,
                                 "Pos": pos,
                                 "Ref": ref,
                                 "Alt": alt,
                                 "RAW": y_pred_prob})
if args.sort:
    df.sort_values("RAW", axis=0, inplace=True, ascending=False)
df.to_csv(args.output, index=False, header=not args.no_header)
