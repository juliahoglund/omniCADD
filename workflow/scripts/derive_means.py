#!/usr/bin/env python3
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
# ...existing code...
