#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 12.03.2018

This files calls up the gerpelem.pl and getgerp.pl and prints out all scores for the requested chromosome. 
The script returns the position, depending on the start position, 
to avoid to add the start position to each position in the output, 
all start positions are requested immediately per chromosome.


:Edited by: Seyan Hu
:Date: 14-11-2022
:Extension and modification: Julia Höglund
:Date: 2023-08-24
:Usage: python <script.py> -o <Outputfile location> -f <Reference chromosome directory> -s <Species of interest>

"""


# Import dependencies. 
import os,sys
from optparse import OptionParser
import pysam


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-o", "--out", dest="out", help="path of the outfile (that will be created)",default="./")
parser.add_option("-f", "--file", dest="file", help="path to reference species chromosome files",default="genome/")
parser.add_option("-s", "--species", dest="species", help="Species of reference",default="Sus scrofa")

(options, args) = parser.parse_args()


# Create list of the reference chromosomes. 
chr_file_list = []
for fn in os.listdir(options.file):
	if fn.endswith('.fa'): 
		chr_file_list.append(fn)

chr_file_list = sorted(chr_file_list)

# Iterate over list.
for chr_ref in chr_file_list:
	
	# Determine chr.
	chr_num = chr_ref.split('_')[-1].replace('.fa', '')

	# Open reference file with pysam.
	ref_fasta = pysam.Fastafile(options.file + chr_ref)
	
	# Determine start and end of chr. 
	# Also checks for the chromosome identifier.
	open_ref = open(options.file + chr_ref, 'r')
	chr_ident = ''
	for line in open_ref:
		line = open_ref.readline().strip()
		# line = line.split(' ')[0]
		line = line.split('>')[1]
		chr_ident = line
		break
	start = 0
	end = ref_fasta.get_reference_length(chr_ident)
	
	
	# Prepartions for GerpElement
	print('Working on GerpElement for Chr: ' + chr_num)
	outfile = open(options.out + chr_num + '_GerpElement.bed','w')
	outfile.write("Chr\tStart\tEnd\tGerpElem\tGerpElemPVal\n")
	outfile.close()
	
	# Calling gerp elements and appending outfile
	os.system("perl workflow/step_5_annotate_variants/scripts/getelem.pl --chr %s --start %s --end %s --sp %s | tail -n +2 >> %s" %(chr_num, str(start), str(end), options.species, options.out + chr_num + '_GerpElement.bed'))
	
	# Prepartions for GerpScore
	print('Working on GerpScore for Chr: ' + chr_num)
	outfile = open(chr_num + '_Gerp_conservation_scores.tsv','w')
	outfile.write("Chr\tStart\tObserved\tExpected\tDifference\n")
	outfile.close()
	
	# Because requesting all data for an entire chromosome does not work for the largest chromosomes, 
	# it is performed with fixed step chunks of 1 mill base pairs.
	start_inner = start
	end_inner = start + 1000000

	while end_inner<=end:
		os.system("perl workflow/step_5_annotate_variants/scripts/getgerp.pl --chr %s --start %s --end %s --sp %s | tail -n +2 | awk '{print $1\"\t\"$2+%s\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}' >> %s" %(chr_num, str(start_inner), str(end_inner), options.species, str(start_inner-1), chr_num + '_Gerp_conservation_scores.tsv'))
		
		start_inner = end_inner
		end_inner += 1000000
	
# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_GERP.txt', 'x')
indication.close()
os.rename('./finished_GERP.txt', './output/finished_GERP.txt')