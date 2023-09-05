#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 5-12-2022
:Extension and modification: Julia Höglund
:Usage: python <script.py> -p <Path to formatted phyloP scores> -f <Path to formatted phastCons scores> -n <List of considered chromosomes>

The script splits the formatted file into subfiles per chromosome. 

:Example:
python PHAST_split.py -f ./phastCons_scores.txt -p phyloP_scores.txt -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X

"""

# Import dependencies.
import sys, os
from os import listdir
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-p","--phylo", dest="phylo", help="path to formatted file containing phyloP scores",default="./phyloP_scores.txt")
parser.add_option("-f","--phast", dest="phast", help="path to formatted file containing phastCons scores",default="./phastCons_scores.txt")
parser.add_option("-n","--chromosome", dest="chromosome", help="list of considered chromosomes",default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X")
parser.add_option("-c","--clean", dest="clean", help="remove original genomewide files? [yes/no]; default: no",default="no")


(options, args) = parser.parse_args()


# Iterate over the chromosomes.  
for chr_num in options.chr.split(','):
	
	# Create output.
	if 'phastCons' in options.phast:
		outp = open(chr_num + '_phastCons.tsv', 'w')
	elif 'phyloP' in options.phylo:
		outp = open(chr_num + '_phyloP.tsv', 'w')
	outp.write('#chr_num\tstart_position\tscore\n')
	
	# Open input.
	open_inp = open(options.file, 'r')
	
	# Iterate over lines in opened file. 
	for line in open_inp:
		split_line = line.split('\t')
		chr_number = split_line[0]
		if chr_number == chr_num:
			outp.write(line)