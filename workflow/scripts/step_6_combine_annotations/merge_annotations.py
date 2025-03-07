#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
script to merge annotations to one annotation file per chromosome.
"""

import pandas as pd
from argparse import ArgumentParser
import logging
from typing import Any

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments() -> Any:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--vep",
        help="Processed file (chromosome wide) with vep annotations", 
        type=str, 
        required=True)
    parser.add_argument("-b", "--bed",
        help="Bed file (chromosome wide) with processed phast and gerp annotations",
        type=str, 
        required=True)
    parser.add_argument("-o", "--outfile",
        help="Name of outfile with combined annotations",
        type=str, 
        required=True)
    return parser.parse_args()

def main() -> None:
    args = parse_arguments()

    try:
        # combine vep and evolutionary constraint
        logging.info("Reading VEP file...")
        vepfile = pd.read_csv(args.vep, sep="\t", low_memory=False)
        
        logging.info("Reading BED file...")
        bedfile = pd.read_csv(args.bed, sep=" ", low_memory=False)
        
        # 1. remove some unwanted columns in bed file, 
        #    like chrom end maybe more
        logging.info("Processing BED file...")
        bedfile = bedfile.drop(columns=['chr', 'end'])
        bedfile = bedfile.rename(columns={"start": "Pos"})
        
        # 2. left outer join with vep versus bed
        logging.info("Merging files...")
        left_merged = pd.merge(vepfile, bedfile, how="left", on=["Pos"])
        
        # 3. write to file
        logging.info("Writing output file...")
        left_merged.to_csv(args.outfile, index=False, sep="\t")
        logging.info("Merge completed successfully.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
