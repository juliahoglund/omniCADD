#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 16-10-2023
:Usage: see derive_means.py --help

Script to derive means from all simulated variants to use for mean
imputation in the dataset. Originally this was part of the data_preparation
script but for optimal parallelization this was split into this script.

Only columns for which the mean is needed, as defined in the annotation
processing configuration, will have be loaded and have their mean stored in
the imputation dictionary.
"""

# Import dependencies
from argparse import ArgumentParser
import pandas
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Create and process 
parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
                    help="Fully annotated input file(s), at least 1 should "
                         "be provided. Additional files will be merged.",
                    type=str, nargs="+", default="chr1_annotations.tsv") # required=True, 
parser.add_argument("-p", "--processing-config",
                    help="Configuration tsv file indicating how the dataset "
                         "should be processed",
                    type=str, default="annot_processing_config.tsv") # required=True, 
parser.add_argument("-o", "--output",
                    help="Dictionary file to write imputation values to.",
                    type=str, default="impute_dict.txt")

args = parser.parse_args()


def load_tsv_configuration(file: str) -> dict[str, dict[str, str]]:
    """
    Loads configuration from tsv table file.
    The first non # or empty line is read as column names

    :param file: str, filename to load configuration from
    :return: dict key is entry label and value is dict,
    with key column label and the value of that column
    """
    try:
        with open(file, "r") as file_h:
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
                samples[parts[0]] = {key: value for key, value in zip(elements, parts[1:])}
    except FileNotFoundError:
        logging.error(f"Configuration file {file} not found.")
        raise
    return samples


def load_mean_cols(infiles: list[str]) -> pandas.DataFrame:
    """
    For each file in the input load pandas dataframe from csv and save
    required columns, as defined in the configuration file
    :param infiles: list of str, input csv files to load
    :return: Pandas dataframe, all files merged and filtered for required cols
    """
    dtypes = dict(
        [(key, value["type"]) for key, value in CONFIGURATION.items()])
    parts = []
    for file in infiles:
        try:
            df = pandas.read_csv(file,
                                 sep='\t',
                                 na_values=['-'],
                                 dtype=dtypes)[MEAN_COLS]
            parts.append(df)
        except FileNotFoundError:
            logging.error(f"Input file {file} not found.")
            raise
        except Exception as e:
            logging.error(f"Error loading file {file}: {e}")
            raise
    return pandas.concat(parts, axis=0)


try:
    # Load configuration from config files
    CONFIGURATION = load_tsv_configuration(args.processing_config)
    MEAN_COLS = [key for key, value in CONFIGURATION.items()
                 if value["isMetadata"] == "False" and value["impute"] == "Mean"]

    # Load data
    myData = load_mean_cols(args.input)

    logging.info("Data from which means are derived:")
    logging.info(myData.describe())

    # Deriving means and writing them to file
    means = {}
    for label in MEAN_COLS:
        means[label] = myData[label].mean()
    with open(args.output, "w") as f:
        f.write(str(means))
    logging.info(f"Means written to {args.output}")
except Exception as e:
    logging.error(f"An error occurred: {e}")