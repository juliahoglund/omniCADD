#!/usr/bin/env python3
# -*- coding: ASCII -*-
"""
:Author: Christian Gross
:Contact: cgross@tudelft.nl
:Date: 01-08-2018

This script takes the vcf file with derived variants and filters out those
who lies on consecutive positions, as those might be interpreted as indels,
i.e. longer than SNPs

:Edited by: Seyan Hu
:Date: 20-10-2022
:Edited by: Job van Schipstal
:Date: 22-9-2023
:Usage: see filter_for_singletons.py --help
Reformatted code in accordance with PEP guidelines and switched to argparse.
"""

# Import dependencies
from argparse import ArgumentParser
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = ArgumentParser(description=__doc__)
parser.add_argument('-i', '--input',
                    help='Input vcf file to filter for SNPs',
                    required=True)
parser.add_argument('--snps',
                    help='Output file path for separated SNPs',
                    required=True)
parser.add_argument('--series',
                    help='Output file for adjacent SNPs',
                    required=True)

# ...existing code...
