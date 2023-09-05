#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 16-11-2022
:Extension and modification: Julia Höglund
:Date 2023-08-21
:Usage: python <script.py> -v <VEP output path> -r <reference chromosome path> -s <original vcf file before VEP annotation path> -g <path to grantham matrix> -t path to VEP_process script>

Wrapper for processing the VEP output on the derived and simulated variants. 

:Example:
python wrapper_VEP_process.py -s output/indexed/ -r genome/ -v output/VEP/ -g data/grantham_matrix/grantham_matrix_formatted_correct.tsv


DOUBLECHECK HARDCODED NAMING OF FILES
"""

# Import dependencies.
import sys,os
from optparse import OptionParser

# OptionParser for input. 
parser = OptionParser()
parser.add_option("-v","--vep", dest = "vep", help = "path to VEP annotation output", default="output/VEP")
parser.add_option("-r","--ref", dest = "ref", help = "path and prefix to reference chr", default = 'genome/Sus_scrofa_ref_')
parser.add_option("-s","--vcf_source", dest="vcf", help="path to bgzipped & tabix index vcf source file (original variant files with simulated and derived variants before annotation)",default="output/indexed/") 
parser.add_option("-g","--grantham", dest = "grantham", help = "Path to Grantham score annotation file", default = 'step_5_annotate_variants/data/grantham_matrix/grantham_matrix_formatted_correct.tsv')
parser.add_option("-t", "--generate", dest="generate", help="path to where the 'VEP_process.py' script is located", default='scripts/')
parser.add_option("-c", "--clean", dest="clean", help="remove previous VEP annotation files? (yes/no; default: no)", default='no')

(options, args) = parser.parse_args()

# Sorted list of files in directory
VEP_list = []
for file_n in os.listdir(options.vep):
	if file_n.endswith('vcf'):
		VEP_list.append(file_n)
VEP_list = sorted(VEP_list)

# Iterate over sorted list.
print('Start processing VEP!')
for fn in VEP_list:
	print('Working on: ' + fn)

	# Check if input is a derived variant. 
	if 'derived' in fn:
		
		# Determine chr number. 
		chr_num_d = fn.split('chr')[-1]
		chr_num_d = chr_num_d.replace('_VEP-annotated.vcf', '')
		
		# Perform commandline.
		os.system('python ' + options.generate + 'VEP_process.py -v ' + options.vep + fn + ' -r ' + options.ref + chr_num_d + '.fa' + ' -s ' + options.vcf + 'derived_variants_chr' + chr_num_d + '_case_upper.vcf.gz' + ' -g ' + options.grantham)
	
	# Check if input is a simulated variant. 
	if 'simulated' in fn:
		
		# Determine chr number. 
		chr_num_d = fn.split('chr')[-1]
		chr_num_d = chr_num_d.replace('_VEP-annotated.vcf', '')
		
		# Perform commandline.
		os.system('python ' + options.generate + 'VEP_process.py -v ' + options.vep + fn + ' -r ' + options.ref + chr_num_d + '.fa' + ' -s ' + options.vcf + 'simSNPs_' + chr_num_d + '.vcf.gz' + ' -g ' + options.grantham)  
		
if options.clean=='yes':
	os.system('rm -r ' + str(options.vep))

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_VEP_processing.txt', 'x')
indication.close()
os.rename('./finished_VEP_processing.txt', './output/finished_VEP_processing.txt')