#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
script to merge annotations to one annotation file per chromosome.
"""

import pandas as pd
from argparse import ArgumentParser

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

args = parser.parse_args()

try:
    # Read VEP file
    try:
        vepfile = pd.read_csv(args.vep, sep="\t", low_memory=False, error_bad_lines=False)
    except pd.errors.ParserError as e:
        raise ValueError(f"Error reading VEP file: {e}")
    
    # Ensure required columns exist in VEP file
    required_vep_columns = {'#Chrom', 'Pos'}
    if not required_vep_columns.issubset(vepfile.columns):
        raise ValueError(f"The VEP file is missing one or more required columns: {required_vep_columns}")
    
    # Rename VEP columns to match BED columns
    vepfile = vepfile.rename(columns={"#Chrom": "chr", "Pos": "start"})
    
    # Read BED file
    try:
        bedfile = pd.read_csv(args.bed, sep=" ", low_memory=False, error_bad_lines=False)
    except pd.errors.ParserError as e:
        raise ValueError(f"Error reading BED file: {e}")
    
    # Ensure required columns exist in BED file
    required_bed_columns = {'chr', 'start', 'end'}
    if not required_bed_columns.issubset(bedfile.columns):
        raise ValueError(f"The BED file is missing one or more required columns: {required_bed_columns}")
    
    # Strip 'chr' prefix from BED 'chr' column
    bedfile['chr'] = bedfile['chr'].str.lstrip('chr')
    
    # Drop unwanted columns in BED file
    bedfile = bedfile.drop(columns=['end'])
    
    # Ensure 'chr' and 'start' column data types match between VEP and BED files
    if vepfile['chr'].dtype != bedfile['chr'].dtype:
        bedfile['chr'] = bedfile['chr'].astype(vepfile['chr'].dtype)
    if vepfile['start'].dtype != bedfile['start'].dtype:
        bedfile['start'] = bedfile['start'].astype(vepfile['start'].dtype)
    
    # Perform left outer join
    left_merged = pd.merge(vepfile, bedfile, how="left", on=["chr", "start"])
    
    # Save output
    left_merged.to_csv(args.outfile, index=False, sep="\t")
    print(f"Successfully merged files and saved to {args.outfile}")

except Exception as e:
    print(f"Error: {e}")
    exit(1)



