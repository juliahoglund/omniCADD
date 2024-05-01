#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 7-12-2023
:Usage: import the helper functions

Provides helper functions for the other scripts.
"""

# Import dependencies
import sys
import os
import pickle
from typing import Union

import pandas
from joblib import Parallel, delayed
from scipy.sparse import load_npz, vstack, save_npz

def filter_matrix(desired, columns, data_matrix):
    """
    Returns desired columns from a scipy data matrix.
    :param desired: list of str, desired columns
    :param columns: list of str, column names of data_matrix
    :param data_matrix: Scipy sparse matrix  to extract columns from
    :return: Pandas dataframe with desired columns
    """
    id_cols = {}
    for col in desired:
        if col not in columns:
            sys.exit(f"Missing column {col} in scipy sparse matrix, "
                     f"but it was requested")
        id_cols[columns.index(col)] = col
    id_cols = dict(sorted(id_cols.items()))
    return data_matrix[:, list(id_cols.keys())], list(id_cols.values())

def read_columns(file: str) -> list:
    """
    Read column csv, return list of columns
    :param file: str, path+file to read columns from
    :return: list of str, columns
    """
    columns = []
    with open(file, "r") as col_f:
        for line in col_f:
            line = line.strip()
            if not line:
                continue
            columns.extend([col.strip() for col in line.split(",")])
    if len(columns) == 0:
        sys.exit(f"Read columns file {columns} and found no columns!")
    return columns

def load_data(file: str, cols: list, desired: Union[list, None], log: bool) \
        -> tuple:
    """
    Loads data in npz format to scipy sparse matrix and metadata in csv
    format to a pandas dataframe, filter the matrix if needed.
    :param file: Iterable of str, files to load the data from.
    :param cols: list of str, column names for the scipy csr
    :param desired: list of str, column names to load (def None -> load all)
    :param log: bool (def True), write progress to stderr
    :return: tuple, scipy csr with desired columns, metadata and selected cols
    """
    selected_cols = None
    if log:
        sys.stderr.write(f"## Loading data from "
                         f"{os.sep.join(file.split(os.sep)[-2:])}\n")

        # Load metadata via pandas df and the actual data via scipy sparse.
    meta_data = pandas.read_csv(file + ".meta.csv.gz", sep=',',
                                na_values=['-'], header=0)
    npz = load_npz(file)
    if desired:
        npz, selected_cols = filter_matrix(desired, cols, npz)
    return npz, meta_data, selected_cols



def load_npz_with_meta(files: list, desired: Union[list, None] = None,
                       log: bool = True, n_jobs: int = 3) -> tuple:
    """
    Loads data in npz format to scipy sparse matrix and metadata in csv
    format to a pandas dataframe, column list is also returned and it is
    verified that this is the same for each file.
    :param files: Iterable of str, files to load the data from.
    :param desired: list of str, column names to load (def None -> load all)
    :param log: bool (def True), write progress to stderr
    :param n_jobs: int (def 3), number of threads to use for reading data
    :return: tuple, scipy sparse matrix, pandas dataframe, list of column names
    """
    cols = None
    for file in files:
        # Load cols and verify they match for the different files
        new_cols = read_columns(file + ".columns.csv")
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

    # Merge if needed
    if len(files) > 1:
        blocks = Parallel(n_jobs=n_jobs)(
            delayed(load_data)(file, cols, desired, log) for file in files
        )
        data = [block[0] for block in blocks]
        meta = [block[1] for block in blocks]
        selected_cols = blocks[0][2]  # Are all the same
        del blocks
        if log:
            sys.stderr.write("## Merging data parts\n")
        # CSR is row based, stacking it by row should be efficient, logistic
        # regression from scikit-learn also expects a csr format sparse matrix
        # so we want it to remain of that type.
        data = vstack(data, format="csr")
        meta = pandas.concat(meta, axis=0)
    else:
        data, meta, selected_cols = load_data(files[0], cols, desired, log)

    if desired:
        cols = selected_cols
    return data, meta, cols


def load_dataset(files, desired_cols, jobs=1):
    """
    Wrapper to load dataset, filter for desired column and take y from metadata
    :param files: list of str, files to read data from
    :param desired_cols: list, desired columns, None for all
    :param jobs: int (def 0 -> len(files)), number of concurrent jobs to run
    :return: tuple, data csr, np array with y values, list of columns for csr
    """
    if jobs < 1:
        jobs = len(files)
    if desired_cols and desired_cols != "All":
        cols = read_columns(desired_cols)
    else:
        cols = None
    data_matrix, metadata, columns = load_npz_with_meta(files,
                                                        cols,
                                                        n_jobs=jobs)
    if "y" in metadata.columns:
        y_values = metadata[["y"]].to_numpy(dtype="float").flatten()
    else:
        y_values = None
    return data_matrix, y_values, columns


def save_npz_with_meta(file, data_m, meta, cols, log=True) -> None:
    """
    Writes sparse data matrix to npz file,
    with metadata and columns in additional csv files.
    :param file: base filename to write to, expected to end in npz
    :param data_m: data matrix, scipy sparse matrix
    :param meta: pandas dataframe, metadata not included in model data matrix
    :param cols: list of str, columns of data matrix
    :param log: bool (def True), write progress to stderr
    :return: None, written to files
    """
    if log:
        sys.stderr.write(f"## Writing data to "
                         f"{os.sep.join(file.split(os.sep)[-2:])}\n")
    save_npz(file, data_m, compressed=True)
    meta.to_csv(file + ".meta.csv.gz", index=False, na_rep="-")
    # Write column names to file, since scipy sparse doesn't support them
    with open(file + ".columns.csv", "w") as f:
        f.write(",".join(cols))


def get_pickle(file: str, clazz):
    """
    Loads an object from file and checks if it is an instance of clazz
    :param file: file to read object from
    :param clazz: clazz that is expected to be returned
    :return: instance of clazz, loaded from file
    """
    with open(file, "rb") as file_h:
        ins = pickle.load(file_h)
    if not isinstance(ins, clazz):
        sys.exit(f"Error loading instance from pickle file, "
                 f"expected an instance of {clazz} but got {type(ins)}")
    return ins