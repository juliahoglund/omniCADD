import argparse
import glob
import os
import pandas as pd
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Combine constraint annotations for a chromosome.")
    parser.add_argument("-c", "--chromosome", type=str, required=True, help="chromosome number")
    parser.add_argument("-f", "--phast_folder", type=str, default="results/annotation/phast/phastCons", help="path to phastCons scores")
    parser.add_argument("-g", "--phylo_folder", type=str, default="results/annotation/phast/phyloP", help="path to phyloP scores")
    parser.add_argument("-i", "--gerp_folder", type=str, default="results/annotation/gerp", help="path to gerp scores")
    parser.add_argument("-o", "--output", type=str, default=None, help="output file name (default: constraint.<chromosome>.bed)")
    return parser.parse_args()

def load_and_concat(folder, chrom, suffix):
# ...existing code...
