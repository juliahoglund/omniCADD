#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 5-12-2022
:Extension and modification: Julia Höglund
:Usage: python <script.py> -f <Path to formatted file containing phastCons/phyloP scores> -c <List of considered chromosomes>

The script splits the formatted file into subfiles per chromosome. 

:Example:
python split_pC_pP_scores.py -f ./phastCons_scores.txt -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X

"""

# Import dependencies.
import sys, os
from os import listdir
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-e","--elem", dest="elem", help="Path to formatted file containing GERP elements scores",default="annotations/GERP_elem.bed")
parser.add_option("-g","--gerp", dest="cons", help="Path to formatted file containing GERP conservation scores",default="annotations/GERP_cons.bed")
parser.add_option("-c","--chr", dest="chr", help="List of considered chromosomes",default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X")
parser.add_option("-t","--tidy", dest="clean", help="Remove big files after reformatting? [yes/no]; default: no", default="no")


(options, args) = parser.parse_args()


# Iterate over the chromosomes.  
for chr_num in options.chr.split(','):
	
	# Create output.
	outp_elem = open(chr_num + '_GERPelem.tsv', 'w')
	outp_elem.write('#chr_num\tstart_position\tscore\n')

	outp_cons = open(chr_num + '_GERPcons.tsv', 'w')
	outp_cons.write('#chr_num\tstart_position\tscore\n')
	
	# Open input.
	open_elem = open(options.elem, 'r')
	
	# Iterate over lines in opened file. 
	for line in open_elem:
		split_line = line.split('\t')
		chr_number = split_line[0]
		if chr_number == chr_num:
			outp_elem.write(line)

	os.system('gzip ' + chr_num + '_GERPelem.tsv')

		# Open input.
	open_cons= open(options.cons, 'r')
	
	# Iterate over lines in opened file. 
	for line in open_cons:
		split_line = line.split('\t')
		chr_number = split_line[0]
		if chr_number == chr_num:
			outp_cons.write(line)

	os.system('gzip ' + chr_num + '_GERPcons.tsv')

if options.clean=='yes':
	os.system('rm ' + str(options.gerp) + ' ' + str(options.elem))
