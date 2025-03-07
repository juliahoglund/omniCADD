#!/usr/bin/env python
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

# Create and process cli
parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="annotated infile containing annotated variants to process", 
                    type=str, 
                    required=True)

parser.add_argument("-c", "--csv",
                    help="output file, processed dataset in csv format (default None, not written)",
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
                    help="configuration tsv file indicating how to dataset should be processed",
                    type=str, 
                    required=True)

parser.add_argument("--interaction-config",
                 help="Configuration tsv file indicating which "
                      "interaction terms should be generated in the dataset",
                 type=str, 
                 required=True)

parser.add_argument("--imputation-dict",
                    help="dictionary file to read/write imputation values to.",
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


######## functions
##################

def load_tsv_configuration(file: str) -> dict:
    """
    Loads configuration from tsv table file.
    The first non # or empty line is read as column names

    :param file: str, filename to load configuration from
    :return: dict key is entry label and value is dict,
    with key column label and the value of that column
    """
    elements = None
    samples = {}
    with open(file, "r") as file_h:
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

def class_encoded_check(classlabel: str, selection: list, data: pandas.DataFrame) -> pandas.DataFrame:
    """
    Checks whether one-hot encoding resulted in all variants, since missing
    columns would pose an issue during when running the model.
    :param classlabel: str, label of original annotation, used as prefix
    :param selection: Iterable, columns to check existence of
    :param data: dataframe to check encoding in
    :return: dataframe that has been encoded
    """
    for clazz in selection:
        col = f"{classlabel}_{clazz}"
        if col not in data:
            data[col] = np.zeros(data.shape[0], dtype=float)
    return data

def is_float(pot_float: str) -> bool:
    """
    Tries to convert str to float to see if a string is a float.
    :param pot_float: str to check for float
    :return: Bool, can be converted to float or not.
    """
    try:
        float(pot_float)
        return True
    except ValueError:
        return False

######## main script
####################

if __name__ == "__main__":
    try:
        # Load configuration from config files
        CONFIGURATION = load_tsv_configuration(args.processing_config)
        ANNOTATIONS = {key: value for key, value in CONFIGURATION.items() if value["isMetadata"] == "False"}
        INTERACTIONS = load_tsv_configuration(args.interaction_config)

        # Build dtype list for opening the dataset:
        dtypes = {key: value["type"] for key, value in CONFIGURATION.items()}

        # Open DF in pandas.
        df = pandas.read_csv(args.input, sep='\t', na_values=['-'], dtype=dtypes)

        #### CHECK but julia does it
        # Rename columns starting with # and switch ref/alt nt/aa for derived variants
        renames = {'Alt': 'Ref', 'Ref': 'Alt', 'nAA': 'oAA', 'oAA': 'nAA'} \
            if args.derived else {}
        for column in df.columns.values:
            if column.startswith("#"):
                renames[column] = column.replace("#", "")
        if len(renames) > 0:
            df.rename(index=str, columns=renames, inplace=True)

        # Verify config and real data columns match
        difference = set(df.columns.values).difference(set(CONFIGURATION.keys()))
        for dif in difference:
            text = F"Mismatch for column {dif}, only occurs in: "
            text += "Configuration" if dif in list(CONFIGURATION.keys()) else "Data"
            sys.stderr.write(text + "\n")
        if difference:
            sys.exit("Exiting since actual data columns and configuration file are mismatched")

        # Loading imputation dict
        means = {}
        try:
            with open(args.imputation_dict, "r") as f:
                means = ast.literal_eval(f.read())
        except IOError:
            sys.exit(
                "ERROR: The imputation_dict.txt file which contains the mean "
                "imputed values from the simulated data is missing. This "
                "information is necessary to continue.")


        # Process the data following the configuration file
        for key, value in ANNOTATIONS.items():
            key = str(key)

            # Create indicator columns before imputation.
            if value["missingIndicator"] == "True":
                df['IND_' + key] = df[key].isnull().astype(float)

            # Impute missing values
            impute = value["impute"]
            if is_float(impute):
                print(f"imputing fixed float value {key}: {float(impute):.2f}")
                if value["type"] == "bool":
                    df[key] = df[key].astype(float).fillna(float(impute))
                else:
                    df[key] = df[key].fillna(float(impute))
            elif impute == "Mean":
                print(f"imputing Mean {key}, {means[key]}")
                try:
                    df[key] = df[key].fillna(means[key])
                except KeyError:
                    sys.exit(f"{args.imputation_dict} was loaded correctly "
                             f"but does not contain all Keys. {key} is missing")
            elif impute != "None":
                sys.stderr.write(f"INFO: using string for imputation, is this correct?"
                                 f"\nKey: {key}, impute: {impute}\n")
                df[key] = df[key].fillna(impute)

            # Remap values in object
            encode = True
            if "remap" in value and not value["remap"].strip() == "":
                print(f"remapping: {key}")
                remap = {}
                pairs = value["remap"].split(",")
                for pair in pairs:
                    try:
                        old, new = pair.split(":")
                        if old == "True":
                            old = True
                        elif old == "False":
                            old = False
                        elif old == "NaN":
                            old = np.nan
                        if is_float(new):
                            new = float(new)
                            encode = False  # Is being directly remapped to float
                        remap[old] = new
                    except ValueError:
                        sys.exit(f"remapping expected to come in pairs old:new\n"
                                 f"got {pair} from annotation {key}")
                df[key] = df[key].map(remap)

            # One hot encoding for object type annotations via pandas.get_dummies
            if value["type"] == "object" and encode:
                df = pandas.get_dummies(df, prefix=key, dtype="float",
                                        columns=[key], sparse=False)

                # Special treatment for consequence, flip stop_gained and lost
                if args.derived and key == "Consequence":
                    columnsTitles = {'Consequence_SG': 'Consequence_SL',
                                     'Consequence_SL': 'Consequence_SG'}
                    df.rename(index=str, columns=columnsTitles, inplace=True)

                # Verify all categories are created, create those not present
                # May not be the case for e.g. a small whole genome CADD block
                categories = set()
                if "remap" in value and not value["remap"].strip() == "":
                    pairs = value["remap"].split(",")
                    for pair in pairs:
                        categories.add(pair.split(":")[1])
                if "checkCategories" in value and \
                        not value["checkCategories"].strip() == "":
                    categories.update(value["checkCategories"].split(","))
                print(categories)
                df = class_encoded_check(key, categories, df)


        # Sorting dataframe columns, to ensure consistent order. Notably the
        # class_encoded_check step, that adds missing columns for when this specific
        # block did not include a possible values will have a different order
        # depending on which columns are present or missing. Therefore we here sort
        # the dataframe, to get a consistent order. The interaction terms will not
        # be uncertain and having them at the end is easier to read so sorting
        # happens before that step.
        df = df.reindex(sorted(df.columns), axis="columns")

        SPARSE_FLOAT = pandas.SparseDtype("float64", fill_value=0)

        # Create interaction terms
        for name, value in INTERACTIONS.items():
            print(f"Generating interactions for {name}")
            columns = [df]
            for i in value["A_cols"].split(","):
                a_col = df[i].values
                for j in value["B_cols"].split(","):
                    b_col = df[j].values
                    columns.append(pandas.Series(
                        [o * a for o, a in zip(a_col, b_col)], name=i + '_' + j,
                        index=df.index, dtype=SPARSE_FLOAT))
            df = pandas.concat(columns, axis=1)


        # If write to csv, do so
        if args.csv != "None":
            print("Writing to csv")
            df.to_csv(args.csv, index=False, na_rep="-")


        # If sparse output is desired, save metadata that doesn't fit in sparse
        # matrix separately and store both files
        if args.npz != "None":
            print("Saving metadata")
            meta_cols = [key for key, value in CONFIGURATION.items()
                         if value["isMetadata"] == "True"]
            print(meta_cols)
            meta_data = df[meta_cols]
            # Add y column
            if args.y_value is not None:
                meta_data.insert(0, "y", args.y_value)
            meta_data.to_csv(args.npz + ".meta.csv.gz", index=False, na_rep="-")
            del meta_data

            print("Converting full dataframe to sparse")
            df = df.drop(columns=meta_cols)
            df = df.astype(SPARSE_FLOAT)
            print(f"Resulting density: {df.sparse.density}")

            # Write column names to file, since scipy sparse doesn't support them
            with open(args.npz + ".columns.csv", "w") as f:
                f.write(",".join(df.columns.values))

            print("Converting to scipy COO matrix")
            coo_data = df.sparse.to_coo()
            del df
            csr_data = csr_matrix(coo_data)
            del coo_data
            print("Saving as csr matrix in npz format")
            save_npz(args.npz, csr_data, compressed=True)

        print("Script completed successfully.")

    except Exception as e:
        print(f"An error occurred: {e}")