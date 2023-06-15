#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Christian Gross
:Contact: cgross@tudelft.nl
:Date: 01-08-2018

This script takes the vcf file with simulated SNPs and iterates over it, 
at each row identifying if the variant is at a position with a known ancestor or not. 
Depending on that it is splitting the file.

:Edited by: Seyan Hu
:Date: 4-11-2022
:Extension and modification: Julia Höglund
:Date 1-5-2023
:Usage: python <script.py>  -i <name and path of indel file> -s <name and path of snp file> -a <ancestor genome path>

:Example:
python filter_simulated.py -i -s -a
"""
# Import dependencies.
import pysam
import sys,os
from optparse import OptionParser

# OptionParser for inputs. 
parser = OptionParser()
parser.add_option("-i", "--indels", dest="indels", help="name of simulated indel variants.", default='./indels_simVariants.vcf')
parser.add_option("-s", "--snps", dest="snps", help="name of simulated SNP variants.", default='snps_simVariants.vcf')
parser.add_option("-a", "--ancestor", dest="ancestor", help="path to ancestor genome files.", default='./output/extracted_ancestor')
(options, args) = parser.parse_args()

# Checking if the path ends with '/' 
if (not options.ancestor.endswith('/')):
	options.ancestor = options.ancestor+'/'

# Create list of chr ancestor files.
anc_l = []
for fn in os.listdir(options.ancestor):
	if fn.endswith('.fa'):
		anc_l.append(fn)

# Open vcf files and ancestor file. 
indel_input = open(options.indels, 'r')
snp_input = open(options.snps, 'r')

# Create output files
outfile1 = open('./' + options.indels.replace(".vcf", "") + '_filtered.vcf','w')	
outfile2 = open('./' + options.snps.replace(".vcf", "") + '_filtered.vcf','w')	

outfile1.write('##fileformat=VCFv4.1\n')
outfile1.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')

outfile2.write('##fileformat=VCFv4.1\n')
outfile2.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')

# Iterate over simulated list and ancestor list.
for anc_file in anc_l:
	print("checking {}".format(anc_file))
	# Open ancestor file. 
	anc_open = open(options.ancestor + anc_file, 'r')
	anc_fasta = pysam.Fastafile(options.ancestor + anc_file)	
	chr_num = ''
	anc_file_n = ''

	for line in anc_open:
		if line.startswith('>'):
			line = line.split(' ')[0]
			line = line.replace('>', '')
			line_n = line.split('.')[0]
			line_chr = line.split('.')[1]
			anc_file_n += line_n
			chr_num += line_chr
			continue

	# reset line search in input
	indel_input.seek(0) 
	# Iterate over lines in indel vcf.
	for lines in indel_input:
		# If lines starts with '#' or 'y' write to output. 
		if (lines[0]=='#') or (lines[0]=='y'):
			next
		# Else split line for chromosome and position of SNP. 
		# And check if SNP in ancestor position. 
		else:
			Chrom = lines.split('\t')[0]
			Pos = lines.split('\t')[1]
			if chr_num == Chrom:

				ancestor = anc_fasta.fetch(anc_file_n + '.' + Chrom, int(Pos) - 1, int(Pos))
			
				if ancestor in 'ACGT':
					outfile1.write(lines)
				else: 
					next
			else:
				next

	# reset line search in input
	snp_input.seek(0) 
	# Iterate over lines in snp vcf. 
	for lines in snp_input:
		# If lines starts with '#' or 'y' write to output. 
		if (lines[0]=='#') or (lines[0]=='y'):
			next
		# Else split line for chromosome and position of SNP. 
		# And check if SNP in ancestor position. 
		else:
			Chrom = lines.split('\t')[0]
			Pos = lines.split('\t')[1]
			if chr_num == Chrom:

				ancestor = anc_fasta.fetch(anc_file_n + '.' + Chrom, int(Pos) - 1, int(Pos))
			
				if ancestor in 'ACGT':
					outfile2.write(lines)
				else: 
					next
			else:
				next

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_vcf_filtering.txt', 'x')
indication.close()
os.rename('./finished_vcf_filtering.txt', './output/finished_vcf_filtering.txt')
