#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to merge annotations to one annotation file per chromosome.
"""

import pandas as pd
from argparse import ArgumentParser
import logging
import sys
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


def read_bed_at_positions(path: str, positions: pd.Index, chunksize: int = 5_000_000) -> pd.DataFrame:
    """
    Read the chromosome-wide constraint BED in chunks, keeping only rows
    whose position occurs in the VEP file.
    """
    def usecols(col: str) -> bool:
        return col not in ("chr", "end")

    kept = [
        chunk[chunk["start"].isin(positions)]
        for chunk in pd.read_csv(path, sep=" ", usecols=usecols, chunksize=chunksize)
    ]
    if not kept:
        return pd.read_csv(path, sep=" ", usecols=usecols, nrows=0)
    return pd.concat(kept, ignore_index=True)


def main() -> None:
    args = parse_arguments()

    try:
        # combine vep and evolutionary constraint
        logging.info("Reading VEP file...")
        vepfile = pd.read_csv(args.vep, sep="\t", low_memory=False)

        logging.info("Reading BED file (only positions present in the VEP file)...")
        bedfile = read_bed_at_positions(args.bed, pd.Index(vepfile["Pos"].unique()))
        bedfile = bedfile.rename(columns={"start": "Pos"})

        # left outer join with vep versus bed
        logging.info("Merging files...")
        left_merged = pd.merge(vepfile, bedfile, how="left", on=["Pos"])

        # write to file
        logging.info("Writing output file...")
        left_merged.to_csv(args.outfile, index=False, sep="\t")
        logging.info("Merge completed successfully.")
    except Exception:
        logging.exception("An error occurred while merging annotations")
        sys.exit(1)


if __name__ == "__main__":
    main()
