#!/usr/bin/env python
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

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help = "variants, npz file(s), at least one should be provided. "
                         "It is expected that there are .meta.csv.gz and columns.csv files present as well",
                    type = str, 
                    required = True, 
                    nargs = "+")
parser.add_argument("-o", "--output",
                    help = "output fold files, the amount of files provided equals the number of folds, "
                           "would not recommend going lower as 4",
                    type = str, 
                    required = True, 
                    nargs = "+")
parser.add_argument("-n", "--n-jobs",
                    help = "int, number of jobs for training the models (default 0, as many as needed, "
                           "limited to half the total cores of the device)",
                    type = int, 
                    default = 0)
parser.add_argument("-m", "--module",
                    help = "location of helper script (optional), not needed if in the same folder",
                    type = str, 
                    required = False)

args = parser.parse_args()

# Load helper functions.
# Snakemake caching makes this more complicated,
# have to enable importing from another folder
if args.module:
    sys.path.append(Path(args.module).parent)
from data_helper import load_npz_with_meta, save_npz_with_meta


def fold_data(infiles, outfiles) -> None:
    """
    Load data from files, merge and write stratified folds.
    :param infiles: list of str, files to load data from
    :param outfiles: list of str, files to write folds to
    :return: None, written to file
    """
    n_jobs = args.n_jobs if args.n_jobs != 0 else len(args.output)

    # Load dataset in parallel
    data_csr, meta_data, columns = load_npz_with_meta(infiles, n_jobs=n_jobs)
    y = meta_data[["y"]].to_numpy(dtype="float").flatten()

    skf = StratifiedKFold(n_splits = len(outfiles),
                          shuffle = True,
                          random_state = 0)

    # write out folds in parallel, notably faster but higher memory usage
    folds = [
        [data_csr[test_index, :], meta_data.iloc[test_index, :]]
        for train_index, test_index in skf.split(data_csr, y)
    ]

    del data_csr, meta_data
    Parallel(n_jobs=n_jobs)(
        delayed(save_npz_with_meta)(
            outfiles[i], data, meta_data, columns
        ) for i, (data, meta_data) in enumerate(folds)
    )


if __name__ == '__main__':
    fold_data(args.input, args.output)
