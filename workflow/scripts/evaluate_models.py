#!/usr/bin/env python3
# -*- coding: ASCII -*-

"""
:Author: Julia Höglund
:Date: 2024
:Usage: evaluate_models.py --help

Evaluates trained models by parsing their statistics files.
Compares different hyperparameter combinations (C and max_iter)
and identifies the best performing model based on ROC-AUC score.
Outputs a summary table and saves the best parameters to JSON.
"""

import sys
import json
import logging
import re
from argparse import ArgumentParser
from pathlib import Path
import pandas as pd

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = ArgumentParser(description=__doc__)
parser.add_argument("--stats",
                    help="Stats files from trained models (.stats.txt)",
                    type=str,
                    required=True,
                    nargs="+")

parser.add_argument("--summary",
                    help="Output file for summary table (TSV format)",
                    type=str,
                    required=True)

parser.add_argument("--best-params",
                    help="Output file for best parameters (JSON format)",
                    type=str,
# ...existing code...
