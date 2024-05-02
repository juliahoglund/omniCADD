#!/usr/bin/env python
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
from argparse import ArgumentParser
from pathlib import Path

import pandas
from joblib import Parallel, delayed

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, confusion_matrix, \
    classification_report

parser = ArgumentParser(description=__doc__)
parser.add_argument("--train",
                    help = "Training dataset, npz file(s), at least one should be provided. It is expected that there are "
                           ".meta.csv.gz and columns.csv files present as well",
                    type = str, 
                    required = True, 
                    nargs = "+")

parser.add_argument("--test",
                    help = "Test dataset, npz file(s), It is expected that "
                           "there are .meta.csv.gz and columns.csv files "
                           "present as well. If none are provided, no testing "
                           "will be done (default is no testing)",
                    type = str, 
                    default = [], 
                    nargs = "*")

parser.add_argument("--columns",
                    help = "file (default None, use all columns), "
                           "only train a model based on the "
                           "columns in this file, csv format",
                    type = str, 
                    default = None)

parser.add_argument("-c",
                    help =  "float, inverse of regularisation strength (default 1.0)."
                            " Models will be trained for each combination "
                            "of iterations and C if multiple values are present.",
                    type = float, 
                    default = [1.0], 
                    nargs = "+")

parser.add_argument("-n", "--n-jobs",
                    help = "int, number of jobs for training the models "
                           "(default 0, as many as needed, limited to half the total cores of the device)",
                    type = int, 
                    default = 0)

parser.add_argument("-i", "--max-iter",
                    help = "integer, maximum number of iterations (default 100). "
                           "Models will be trained for each combination of iterations and C if multiple values are present.",
                    type = int, 
                    default = [100], 
                    nargs = "+")

parser.add_argument("--file-pattern",
                    help = "File pattern to write model to. pickle format will be appended. "
                           "must inlcude the tags [C] and [ITER] for the regulization strength "
                           "and maximum number of training iterations, respectively. "
                           "(default model_[C]C_[ITER]iter.pickle)",
                    type = str, 
                    default = "model_[C]C_[ITER]iter.pickle")

parser.add_argument("--save-scaler",
                    help = "File to write sklearn StandardScaler used to scale the data to. "
                           "(default None, not written)",
                    type = str, 
                    default = "None")

parser.add_argument("-w", "--save-weights",
                    help = "bool (default False, don't write), "
                           "Write weights to a csv file",
                    action = "store_true")

parser.add_argument("-m", "--module",
                    help = "Location of helper script (optional), "
                           "not needed if in the same folder",
                    type = str, 
                    required = False)

args = parser.parse_args()

# Load helper functions.
# Snakemake caching makes this more complicated,
# have to enable importing from another folder
if args.module:
    sys.path.append(Path(args.module).parent)
from data_helper import load_dataset


def fit_model(c, max_iter, train_x, train_y):
    """
    Fit a logistic regression model on the training dataset.
    :param c: float, inverse of regularisation strength
    :param max_iter: int, model parameter, max training iterations
    :param train_x: matrix of float, x values to train on
    :param train_y: array of float, y values, 'true' y
    :return: trained instance of LogisticRegression model
    """
    reg_model = LogisticRegression(penalty="l2", solver="lbfgs",
                                   C=c, max_iter=max_iter)
    reg_model.fit(train_x, train_y)
    return reg_model


def get_file(c, max_iter):
    """
    Get filename based on the given pattern.
    :param c: float, inverse of regularisation strength
    :param max_iter: int, model parameter, max training iterations
    :return: str, filled in filename
    """
    out_file = args.file_pattern.replace("[C]", str(c))
    return out_file.replace("[ITER]", str(max_iter))


if __name__ == '__main__':
    # Load train dataset
    x_train, y_train, columns = load_dataset(args.train,
                                             args.columns,
                                             jobs=args.n_jobs)

    sys.stderr.write("## Scaling dataset\n")
    # Mean scaling would break sparsity, causing 100's of Gigabytes of
    # memory usage
    scaler = StandardScaler(with_mean=False, with_std=True, copy=False)
    scaler.fit(x_train)
    # print(scaler.mean_)

    x_train = scaler.transform(x_train)
    # Save scaler if asked
    if args.save_scaler != "None":
        with open(args.save_scaler, "wb") as handle:
            pickle.dump(scaler, handle)

    # Create parameter space
    params = []
    for strength in args.c:
        for iter_count in args.max_iter:
            params.append([strength, iter_count])

    sys.stderr.write("## Fitting Model(s)\n")
    n_jobs = args.n_jobs
    if n_jobs == 0:
        # Then try all in parallel, limit to half the cpu
        n_jobs = min(len(params), int(os.cpu_count() / 2))
    if n_jobs == 1:
        # Save performance/memory overhead of parallel, just loop
        models = [fit_model(c, max_iter, x_train, y_train)
                  for c, max_iter in params]
    else:
        # This runs the training in parallel, in a fairly memory efficient
        # manner, since the dataset is put in shared memory.
        models = Parallel(n_jobs=n_jobs)(delayed(fit_model)(
            c, max_iter, x_train, y_train
        ) for c, max_iter in params)

    # Write models to disk
    sys.stderr.write("## Saving model(s) to disk\n")
    for model, (c, max_iter) in zip(models, params):
        with open(get_file(c, max_iter) + ".pickle", "wb") as model_f:
            pickle.dump(model, model_f)

    if args.save_weights:
        for model, (c, max_iter) in zip(models, params):
            print(f"Length coef-cols: {len(model.coef_[0])}-{len(columns)}")
            with open(get_file(c, max_iter) + ".weights.csv", "w") as weight_f:
                weight_f.write("name,weight\n")
                for coef, name in zip(model.coef_[0], columns):
                    weight_f.write(f"{name},{coef}\n")

    # Delete train dataset to free up memory
    del x_train, y_train

    # Predict on test-set if given
    if len(args.test) == 0:
        sys.exit(0)
    x_test, y_test, columns = load_dataset(args.test,
                                           args.columns,
                                           jobs=args.n_jobs)
    x_test = scaler.transform(x_test)

    for model, (c, max_iter) in zip(models, params):
        sys.stderr.write(f"## Saving tests for model c={c}, "
                         f"max_iter={max_iter}\n")
        y_pred = pandas.Series(model.predict(x_test))
        y_pred_prob = pandas.Series(model.predict_proba(x_test)[::, 1])
        result = pandas.concat({"y_true": pandas.Series(y_test),
                                "y_pred": y_pred,
                                "y_pred_prob": y_pred_prob},
                               axis=1)
        # Save predicted values
        result.to_csv(get_file(c, max_iter) + ".pred.csv.gz", index=False,
                      float_format="%.10g")

        # Write basic statistics on test set
        with open(get_file(c, max_iter) + ".stats.txt", "w") as stat_f:
            stat_f.write(f"Stats for c={c}, max_iter={max_iter}\n")
            stat_f.write("Confusion matrix:\n")
            stat_f.write(str(confusion_matrix(y_test, y_pred)))
            stat_f.write("\nStandard scikit-learn classification "
                         "report output:\n")
            stat_f.write(str(classification_report(y_test, y_pred)))
            auc = roc_auc_score(y_test, y_pred_prob)
            stat_f.write(f"\nROC-AUC score: {auc}\n")
