#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Christian Gross
:Contact: cgross@tudelft.nl
:Date: 01-08-2018

:Edited by: Seyan Hu
:Date: 20-10-2022
:Extension and modification: Julia Höglund

:Usage: python <python file> -p <path to vcf file> -i <name of vcf file>

:Example:
python split_vcf.py -p ./ -i simVariants.vcf

"""

# Import dependencies
import os,sys
from optparse import OptionParser

# OptionParser for the previously generated derived variants files (vcf)
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to derived variants.", default='simVariants.vcf') 
parser.add_option("-p", "--path", dest="path", help="path to folder with vcf files.", default= "./") 
(options, args) = parser.parse_args()

# Checking if the path ends with '/' 
if (not options.path.endswith('/')):
	options.path = options.path+'/'

# Read vcf file and create two seperate files for snps and indels.
infile = open(options.path + options.input,'r')
outfile1 = open(options.path + 'snps_' + options.input,'w')
outfile2 = open(options.path + 'indels_' + options.input,'w')

outfile1.write('##fileformat=VCFv4.1\n')
outfile1.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')

outfile2.write('##fileformat=VCFv4.1\n')
outfile2.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')

# Loop through the lines in vcf file
variants = []
for lines in infile:
		
	# Skips header
	if lines.startswith('#'):
		continue
	
	line = lines.split('\t')
	if len(line[3]) > 1 or len(line[4]) > 1:
		outfile2.write(lines)
	else:
		outfile1.write(lines)

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_vcf_splitting.txt', 'x')
indication.close
os.rename('./finished_vcf_splitting.txt', './output/finished_vcf_splitting.txt')


