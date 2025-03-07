#!/usr/bin/env python3
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 5-12-2022
:Usage: python <script.py> 

These fasta files containing repeats should be downloaded from the UCSC Genome Browser database for the species of interest. 
And should be put in the 'data/repeats/' directory. 
UCSC should have masked fasta files (repeats are in lower case) per chromosome for the species of interest. 
These files should be decompressed. 
The script creates per chromosome a output file containing a list of the position of its repeats. 

:Example:
python get_position_repeats.py -r ./

"""

# Import dependencies.
import sys, os
from os import listdir
import argparse
import json

# ArgumentParser for input.
parser = argparse.ArgumentParser(description='Process fasta files to find repeat positions.')
parser.add_argument("-r", "--repeats", dest="repeats", help="Path to repeats", default="./")

args = parser.parse_args()

# Check if the input directory exists
if not os.path.exists(args.repeats):
    print(f"Error: The directory {args.repeats} does not exist.")
    sys.exit(1)

# Define input files.
fn_list = [fn for fn in os.listdir(args.repeats) if fn.endswith('.fa')]

# Check if any .fa files are found
if not fn_list:
    print(f"No .fa files found in the directory {args.repeats}.")
    sys.exit(1)

# Ensure the output directory exists
output_dir = './output/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Iterate over the file list.
for fn in fn_list:
    print('Working on: ' + fn)
    
    # Determine the chr number.
    chr_num = fn.replace('.fa', '')
    
    # Open file.
    try:
        with open(os.path.join(args.repeats, fn), 'r') as open_fn:
            # Skip the first two lines
            open_fn.readline() # empty line
            open_fn.readline() # header
            # Read the entire sequence
            seq_string = open_fn.read().replace('\n', '')
    except IOError as e:
        print(f"Error opening file {fn}: {e}")
        continue

    # Iterate over position in sequence and add positions of repeats to list.
    print('Creating list of repeat positions for: ' + fn)
    repeat_pos = [pos + 1 for pos, nt in enumerate(seq_string) if nt.islower()]
    
    # Create output with json.
    try:
        with open(os.path.join(output_dir, f"repeats_chr{chr_num}.txt"), "w") as fp:
            json.dump(repeat_pos, fp)
    except IOError as e:
        print(f"Error writing file repeats_chr{chr_num}.txt: {e}")
        continue

# Create a txt file indicating that this process is finished (for snakemake)
try:
    with open('finished_get_repeat_position.txt', 'x') as indication:
        pass
    os.rename('finished_get_repeat_position.txt', os.path.join(output_dir, 'finished_get_repeat_position.txt'))
except IOError as e:
    print(f"Error creating finished indication file: {e}")