#!/usr/bin/env python
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
from optparse import OptionParser
import json

# OptionParser for input. 
parser = OptionParser()
parser.add_option("-r","--repeats", dest="repeats", help="Path to repeats", default = "./")

(options, args) = parser.parse_args()


# Define input files.
fn_list = []
for fn in os.listdir(options.repeats):
	if fn.endswith('.fa'):
		fn_list.append(fn)


# Iterate over the file list.
for fn in fn_list:
	print('Working on: ' + fn)
	
	# Determine the chr number.
	chr_num = fn.replace('.fa', '')
	
	# Open file.
	open_fn = open(options.repeats + fn, 'r')
	
	seq_string = open_fn.readline() # empty line
	seq_string = open_fn.readline() # header
	seq_string = open_fn.readline() # sequence


	# Iterate over position in sequence and add positions of repeats to list.
	print('Creating list of repeat positions for: ' + fn)
	repeat_pos = [] 
	for pos, nt in enumerate(seq_string):
		pos += 1
		if nt.islower():
			repeat_pos.append(pos)
	
	# Create output with json. 
	with open("repeats_chr" + str(chr_num) + '.txt', "w") as fp:
		json.dump(repeat_pos, fp)
		
	
# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_get_repeat_position.txt', 'x')
indication.close()
os.rename('./finished_get_repeat_position.txt', './output/finished_get_repeat_position.txt')