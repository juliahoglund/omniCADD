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
parser.add_option("-p","--phylo", dest="phylo", help="path to formatted file containing phyloP scores",default="./phyloP_scores.bed")
parser.add_option("-f","--phast", dest="phast", help="path to formatted file containing phastCons scores",default="./phastCons_scores.bed")
parser.add_option("-n","--chromosomes", dest="chromosomes", help="list of considered chromosomes",default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X")
parser.add_option("-c","--clean", dest="clean", help="remove original genomewide files? [yes/no]; default: no",default="no")


(options, args) = parser.parse_args()

os.system('gunzip ' + str(options.phast))
os.system('gunzip ' + str(options.phylo))

open_pC = open(options.phast.replace('.gz', ''), 'r')
open_pP = open(options.phylo.replace('.gz', ''), 'r')

# Iterate over the chromosomes.  
for chr_num in options.chromosomes.split(','):

	# phastCons
	outp_pC = open('chr' + chr_num + '_phastCons.bed', 'w')
	outp_pC.write('#chr_num start_position end_position score\n')

	# Iterate over lines in opened file. 
	for line in open_pC:
		split_line = line.split('\t')
		chr_number = split_line[0].replace('chr', '')
		if chr_number == chr_num:
			outp_pC.write(line)
	
	os.system('gzip chr' + str(chr_num) + '_phastCons.bed')

	# phyloP
	outp_pP = open('chr' + chr_num + '_phyloP.bed', 'w')
	outp_pP.write('#chr_num\tstart_position\tend_position\tscore\n')
	
	# Iterate over lines in opened file. 
	for line in open_pP:
		split_line = line.split('\t')
		chr_number = split_line[0].replace('chr', '')
		if chr_number == chr_num:
			outp_pP.write(line)

	os.system('gzip chr' + str(chr_num) + '_phyloP.bed')


os.system('gzip ' + str(options.phylo.replace('.gz', '')))
os.system('gzip ' + str(options.phast.replace('.gz', '')))