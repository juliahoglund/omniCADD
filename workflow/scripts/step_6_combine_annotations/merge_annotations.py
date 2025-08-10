#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
script to merge annotations to one annotation file per chromosome.
"""

import pandas as pd
import gzip
import pysam
from argparse import ArgumentParser

def open_file_for_pandas(filename):
    """
    Open a file that can be either regular, gzipped, or bgzipped for pandas reading.
    Returns the appropriate file handle that pandas can use.
    """
    if filename.endswith('.gz'):
        # Try bgzipped first (common for bioinformatics files)
        try:
            # Test if it's bgzipped by trying to open with pysam
            test_file = pysam.TabixFile(filename, mode='r')
            test_file.close()
            # If successful, it's bgzipped - use gzip for pandas (works for both)
            return gzip.open(filename, 'rt')
        except:
            # If pysam fails, try regular gzip
            try:
                return gzip.open(filename, 'rt')
            except:
                # If both fail, treat as regular file
                return open(filename, 'r')
    else:
        # Check if it's a compressed file without .gz extension
        try:
            with open(filename, 'rb') as f:
                magic = f.read(2)
                if magic == b'\x1f\x8b':  # gzip magic number
                    return gzip.open(filename, 'rt')
        except:
            pass
        # Regular file
        return open(filename, 'r')

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v", "--vep",
    help="Processed file (chromosome wide) with snpeff/vep annotations", 
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
    # Read constraint/bed file (from combine_constraint_anno.py) fully into memory (usually smaller)
    try:
        with open_file_for_pandas(args.bed) as bed_file:
            bedfile = pd.read_csv(bed_file, sep="\t", low_memory=False)
    except pd.errors.ParserError as e:
        raise ValueError(f"Error reading constraint/bed file: {e}")

    # Accept both '#Chrom' or 'chr' for chromosome, and 'Pos' for position
    bed_col_map = {}
    for col in bedfile.columns:
        if col.lower() in ['#chrom', 'chr']:
            bed_col_map[col] = '#Chrom'
        if col.lower() == 'pos':
            bed_col_map[col] = 'Pos'
    bedfile = bedfile.rename(columns=bed_col_map)

    required_bed_columns = {'#Chrom', 'Pos'}
    if not required_bed_columns.issubset(bedfile.columns):
        raise ValueError(f"The constraint/bed file is missing one or more required columns: {required_bed_columns}")

    bedfile['#Chrom'] = bedfile['#Chrom'].astype(str).str.replace('chr', '', regex=False).astype(int)
    bedfile['Pos'] = bedfile['Pos'].astype(int)

    # Set index on bedfile for faster join (if memory allows)
    bedfile.set_index(['#Chrom', 'Pos'], inplace=True)

    # Process VEP/SnpEff file in chunks for memory efficiency
    chunk_size = 50000  # Reduce chunk size for lower memory usage
    first_chunk = True
    
    # Open the VEP file with compression support
    vep_file = open_file_for_pandas(args.vep)
    
    try:
        for vep_chunk in pd.read_csv(vep_file, sep="\t", low_memory=False, chunksize=chunk_size):
            vep_col_map = {}
            for col in vep_chunk.columns:
                if col.lower() in ['#chrom', 'chrom']:
                    vep_col_map[col] = '#Chrom'
                if col.lower() == 'pos':
                    vep_col_map[col] = 'Pos'
            vep_chunk = vep_chunk.rename(columns=vep_col_map)

            required_vep_columns = {'#Chrom', 'Pos'}
            if not required_vep_columns.issubset(vep_chunk.columns):
                raise ValueError(f"The SnpEff/VEP file is missing one or more required columns: {required_vep_columns}")

            vep_chunk['#Chrom'] = vep_chunk['#Chrom'].astype(str).str.replace('chr', '', regex=False).astype(int)
            vep_chunk['Pos'] = vep_chunk['Pos'].astype(int)

            # Set index for join
            vep_chunk.set_index(['#Chrom', 'Pos'], inplace=True)

            # Join using index for efficiency
            merged = vep_chunk.join(bedfile, how="left", rsuffix='_bed').reset_index()

            # Write header only for the first chunk
            merged.to_csv(args.outfile, index=False, sep="\t", mode='w' if first_chunk else 'a', header=first_chunk)
            first_chunk = False
    finally:
        vep_file.close()

    print(f"Successfully merged files and saved to {args.outfile}")

except Exception as e:
    print(f"Error: {e}")
    exit(1)



