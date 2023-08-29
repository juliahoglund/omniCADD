#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 16-11-2022
:Extension and modification: Julia Höglund
:Date: 2023-08-17
:Usage: python <script.py> -v <path to simulated variants> -d <path to derived variants> -s <species of interest>

Performes VEP annotations on the derived and simulated variants. 

:Example: 
python VEP_annotate.py -v output/simulated/ -d output/derived/ -s sus_scrofa
"""

# Import dependencies.
import math
import sys,os
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-v", "--vcf", dest="vcf", help="Name and path to input files", default= "output/genomeWide/")
parser.add_option("-s", "--species", dest="species", help="Species of interest", default= "sus_scrofa")
(options, args) = parser.parse_args()


# Sorted list of files in directory
file_list = os.listdir(options.vcf)
file_list = sorted(file_list)

# Iterate over sorted list.
print('Performing VEP!')
for fn in file_list:
	print('Working on: ' + fn) 
	
	# Determine chr number. 
	chr_num = fn.split('_all_variants')[0]
	chr_num = chr_num.split('_')[-1]

	# unzip vcf
	os.system('gunzip ' + options.vcf + fn)
	fn = fn.replace('.gz', '')

	# Perform commandline.
	# server with cache dir: add --dir $VEP_CACHE remove --offline
	os.system('vep --input_file '+ options.vcf + fn +' --quiet --cache --dir $VEP_CACHE --buffer 1000 --no_stats --species '+ options.species +' --format vcf --regulatory --sift b --per_gene --ccds --domains --numbers --canonical --total_length --force_overwrite --output_file temp.vcf')
	os.system('''cat temp.vcf | awk 'BEGIN{ FS="\t"; OFS="\t"; }{ if ($1 ~ /^#/) { if ($1 ~ /^#Up/) { sub("#","",$1); print "#Chrom","Start","End",$0 } else { print } } else { split($2,a,":"); split(a[2],b,"-"); if (length(b) == 2) { print a[1],b[1],b[2],$0 } else { print a[1],b[1],b[1],$0 } }}' >> allVars_chr'''+ chr_num +'''_VEP-annotated.vcf''')
	
	# zip vcf back
	os.system('gzip ' + options.vcf + fn)
	# zip created output (is bgzip better already here?)
	os.system('gzip allVars_chr' + chr_num + '_VEP-annotated.vcf')


# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_VEP_all_variants.txt', 'x')
indication.close()
os.rename('./finished_VEP_all_variants.txt', './output/finished_VEP_all_variants.txt')

