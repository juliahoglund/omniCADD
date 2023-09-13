#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 12-12-2022
:Extension and modification: Julia Höglund
:Date: 2023-09-07
:Usage: python <script.py> -g <Path to gerpElem output>

Formats the GerpElem file so that it is easier to merge without the use of dicts 
(previously converted to a dict before merging, but uses too much memory). 

:Example:
python GERP_format.py -g annotations/GERP_elem.bed

"""

# Import dependencies. 
import sys, os
from os import listdir
from os.path import isfile, join
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-g","--gerp", dest="gerpE", help="Path to GERP elem file", default="annotations/GERP_elem.bed")
parser.add_option("-t","--tidy", dest="clean", help="Remove wig files after reformatting? [yes/no]; default: no", default="no")


(options, args) = parser.parse_args()
		

# Open file.
open_fn = open(options.gerpE, 'r')
	
	
# Create output.
output = open('formatted_GerpElem.tsv', 'w')
output.write('#Chr\tStart\tGerpElem\n')
	
# Iterate over lines in file.
for line in open_fn:
	split_ln = line.replace('\n', '').split('\t')
	
	chr_num = split_ln[0]
	start_pos = split_ln[1]
	end_pos = split_ln[2]
	GE_score = split_ln[4] # Value, in downloaded, can be 3 if retrieved by API
	for pos in range(int(start_pos), int(end_pos) + 1):
		output.write(str(chr_num) + '\t' + str(pos) + '\t' + str(GE_score) + '\n')

if options.clean=='yes':
	os.system('rm ' + str(options.gerp))
