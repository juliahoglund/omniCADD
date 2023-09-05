#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 12-12-2022
:Usage: python <script.py> -g <Path to gerpElem output>

Formats the GerpElem file so that it is easier to merge without the use of dicts 
(previously converted to a dict before merging, but uses too much memory). 

:Example:
python GERP_format.py -g ./

"""

# Import dependencies.
import sys, os
from os import listdir
from os.path import isfile, join
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-g","--gerpE", dest="gerpE", help="Path to gerpElem output",default="output/")

(options, args) = parser.parse_args()


# List of files. 
gE_list = []
for fn in os.listdir(options.gerpE):
	if fn.endswith('_GerpElement.bed'): # change later only element not conservation?
		gE_list.append(fn)
gE_list = sorted(gE_list)
		

# Iterate over list.
for fn in gE_list:
	
	# Open file.
	open_fn = open(options.gerpE + fn, 'r')
	
	# Determine chr number. 
	chr_num = fn.replace('_GerpElement.bed', '')
	
	# Create output.
	output = open(str(chr_num) + '_formatted_GerpElem.tsv', 'w')
	output.write('#Chr\tStart\tGerpElem\n')
	
	# Iterate over lines in file.
	for line in open_fn:
		if not line.startswith('Chr'):
			split_ln = line.replace('\n', '').split('\t')
			
			start_pos = split_ln[1]
			end_pos = split_ln[2]
			GE_score = split_ln[3] # Value
			for pos in range(int(start_pos), int(end_pos) + 1):
				output.write(str(chr_num) + '\t' + str(pos) + '\t' + str(GE_score) + '\n')