#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 5-12-2022
:Extension and modification: Julia Höglund
:Date: 2023-09-05
:Usage: python <script.py> -w <Path to converted wig file (from bigwig)> -l <path and name of liftover chain file> -t <clean yes/no>

These PhastCons & PhyloP scores should be downloaded from the UCSC Genome Browser database for the species of interest (bigwig format). 
And should be put in the 'data/phastCons/' and 'data/phyloP/' directory. 
The script reformats these scores so that it is easier to parse to the processed VEP output. 
It outputs two files containing the chromosome number, position and score.

:Example:
python PHAST_liftOver.py -w ./ -l annotations/hg38ToSusScr11.over.chain.gz -t no

# https://genome.ucsc.edu/cgi-bin/hgTables
# https://hgdownload.soe.ucsc.edu/downloads.html
# http://hgdownload.cse.ucsc.edu/goldenpath/hg38/
"""

# Import dependencies.
import sys, os
from os import listdir
from optparse import OptionParser
import json


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-p","--phast", dest="phast", help="Path to formatted PHAST files", default="annotations/")
parser.add_option("-l","--liftover", dest="liftover", help="Path (and name of) lift-over file", default="annotations/hg38ToSusScr11.over.chain.gz")
parser.add_option("-t","--tidy", dest="clean", help="Remove wig files after reformatting? [yes/no]; default: no", default="no")

(options, args) = parser.parse_args()

# 1. unzip
# 2 liftover
# 3 zip
# 4. remove if yes 
filenames = [x for x in os.listdir(options.phast) if x.endswith('scores.bed.gz')]

for file in filenames:
	print('working on ' + str(options.phast) + str(file))
	os.system('gunzip ' + str(options.phast) + str(file))
	file = file.replace('.gz', '')
	os.system('liftOver ' + str(options.phast) + str(file) + ' ' + str(options.liftover) + ' ' + str(file.replace('.bed', '_lifted.bed')) + ' unMapped_' + str(file))
	os.system('gzip ' + str(file.replace('.bed', '_lifted.bed')))
	os.system('gzip ' + str(options.phast) + str(file))

if options.clean=='yes':
	os.system('rm ' + str(options.phast) + '*phastCons_scores.bed.gz')
	os.system('rm ' + str(options.phast) + '*phyloP_scores.bed.gz')

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_phast_lifover.txt', 'x')
indication.close()
os.rename('./finished_phast_lifover.txt', './output/finished_phast_lifover.txt')