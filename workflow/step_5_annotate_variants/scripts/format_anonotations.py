#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 5-12-2022
:Usage: python <script.py> -w <Path to converted wig file (from bigwig)>

These PhastCons & PhyloP scores should be downloaded from the UCSC Genome Browser database for the species of interest (bigwig format). 
And should be put in the 'data/phastCons/' and 'data/phyloP/' directory. 
The script reformats these scores so that it is easier to parse to the processed VEP output. 
It outputs two files containing the chromosome number, position and score.

:Example:
python format_pC_pP_wig.py -w ./ -c no

# https://genome.ucsc.edu/cgi-bin/hgTables
# https://hgdownload.soe.ucsc.edu/downloads.html
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/


"""

# Import dependencies.
import sys, os
from os import listdir
from optparse import OptionParser
import json


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-w","--wig", dest="wig", help="Path to converted wig files", default="annotations/")
parser.add_option("-c","--clean", dest="clean", help="Remove wig files after reformatting? [yes/no]; default: no", default="no")


(options, args) = parser.parse_args()


# Defile phastCons, phyloP files. 
for fn in os.listdir(options.wig):
	if fn.startswith('phastCons'):
		pC_fn = fn
	elif fn.startswith('phyloP'):
		pP_fn = fn
		
# Function for creating a dict for the scores in the files.
# Input: Opened phastCons or phyloP wig file with its output name. 
# Output: None and an output file containing the chr number, position and score. 

def scoresToDict(open_f, outn):
	# Create output. 
	outp = open(outn, 'w')
	outp.write('#chr_num\tstart_position\tscore\n')
	
	# Iterate over the lines
	header_sig = 1
	chr_num = 0
	start = 0
	step = 0
	for line in open_f:
		if line.startswith('fixedStep'):
			# Reset stats.
			header_sig = 1
			chr_num = 0
			start = 0
			step = 0
			
			line_split = line.split(' ')
			chr_num = line_split[1].split('=chr')[1]
			start = int(line_split[2].split('=')[1])
			step = int(line_split[3].split('=')[1])
			#print(chr_num, start, step)
			header_sig = 0
			
		elif header_sig == 0:
			score = line.replace('\n', '')
			outp.write(str(chr_num)+'\t'+str(start)+'\t'+str(score)+'\n')
			start += step
			
	return None


# Open files.
open_pC = open(options.wig + pC_fn, 'r')
open_pP = open(options.wig + pP_fn, 'r')


# Perform function 'scoresToDict'.
print('Working on ' + pC_fn)
pC_out_dict = scoresToDict(open_pC, "phastCons_scores.txt")
print('Working on ' + pP_fn)
pP_out_dict = scoresToDict(open_pP, "phyloP_scores.txt")

if options.clean=='yes':
	os.system('rm ' + str(options.wig) + 'phastCons*')
	os.system('rm ' + str(options.wig) + 'phyloP*')


# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_format_wig.txt', 'x')
indication.close()
os.rename('./finished_format_wig.txt', './output/finished_format_wig.txt')