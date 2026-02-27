#!/usr/bin/env python3
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 18-02-2023
:Usage: see --help

Assess column correlation and relevance.
"""

# Import dependencies
from argparse import ArgumentParser
import pandas
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Create and process cli
parser = ArgumentParser(description=__doc__)
parser.add_argument("-d", "--derived",
                    help="annotated infile(s) containing derived variants to process", 
                    nargs="+", 
                    required=True)
parser.add_argument("-s", "--simulated",
                    help="annotated infile(s) containing simulated variants to process", 
                    nargs="+", 
                    required=True)
parser.add_argument("-o", "--outfolder",
                    help="Output folder for correlation tsv files",
                    type=str, 
                    default="None")

args = parser.parse_args()


def get_dataset(infiles: list) -> pandas.DataFrame:
# ...existing code...
