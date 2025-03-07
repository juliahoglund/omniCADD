#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
:Author: Job van Schipstal
:Date: 04-10-2023
:Usage: python <script.py> -i <input.vcf> -o <output.vcf>
 -c <n_existing> -d <n_desired>

Trims the vcf file to the desired number of variants. This is done because
the for training the model we desire the same amount of derived and
ancestral variants, but the amount of simulated variants that will result
from the simulation and subsequent filtering is not exactly known. Trimming
is performed in-place with the use of a mask determining which rows will be
kept and which are discarded.
"""

# Import dependencies.
import sys
import gzip
from argparse import ArgumentParser
import numpy as np
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

parser = ArgumentParser(description=__doc__)
parser.add_argument('-i', '--input',
    help='Input vcf file of variants to trim',
    type=str, 
    required=True)
parser.add_argument('-o', '--output',
    help='output filename, vcf format, can be gzipped',
    type=str,
    required=True)
parser.add_argument('-c', '--current',
    help='Current amount of variants',
    type=str, 
    required=True)
parser.add_argument('-d', '--desired',
    help='Desired amount of variants',
    type=str, 
    required=True)


def trim_vcf(in_h: 'file', out_h: 'file', current: int, desired: int) -> None:
    """
    Trims input vcf file of #current variants to output vcf with #desired
    variants. Filtering is performed in-place with the use of a mask
    determining which rows will be kept and which discarded.

    :param in_h: input vcf file handle.
    :param out_h: output file handle.
    :param current: int, amount of variants in input vcf to trim.
    :param desired: int, desired amount of variants in output vcf.
    :return: None, written to file.
    """
    logging.debug(f"Current number of variants: {current}")
    logging.debug(f"Desired number of variants: {desired}")

    # Get randomised mask, zeros is faster than directly creating boolean mask
    mask = np.zeros(current, dtype=int)
    mask[:desired] = 1
    np.random.shuffle(mask)
    mask = mask.astype(bool)
    logging.debug(f"Mask after shuffling: {mask}")

    obtained = 0
    entry = 0
    for line in in_h:
        if line.startswith("#"):
            out_h.write(line)
            continue
        if mask[entry]:
            out_h.write(line)
            obtained += 1
        entry += 1

    logging.debug(f"Number of variants obtained: {obtained}")

    # Verify the right amount of variants was obtained
    if obtained != desired:
        sys.exit(f"Error: Insufficient variants were output: {obtained} of "
                 f"the desired {desired}, from source file with {current}")


if __name__ == '__main__':
    args = parser.parse_args()
    logging.debug(f"Input file: {args.input}")
    logging.debug(f"Output file: {args.output}")

    # Double conversion enables parsing of scientific notation numbers
    current_n = int(float(args.current))
    desired_n = int(float(args.desired))
    if desired_n > current_n:
        sys.exit(f"Error: Desired number of variants ({desired_n}) exceeds current number of variants ({current_n})."
                 f"\nFile: {args.input}")

    with (gzip.open(args.input, "rt") if args.input.endswith('.gz') else open(args.input, "r")) as vcf_file, \
         (gzip.open(args.output, "wt") if args.output.endswith('.gz') else open(args.output, "w")) as outfile:
        trim_vcf(vcf_file, outfile, current_n, desired_n)

    logging.debug("VCF trimming completed successfully.")
