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
parser.add_option("-v", "--simulated", dest="simulated", help="Name and path of input files", default= "output/simulated/")
parser.add_option("-d", "--derived", dest="derived", help="Name and path of input files", default= "output/derived/")
parser.add_option("-s", "--species", dest="species", help="Species of interest", default= "sus_scrofa")
(options, args) = parser.parse_args()


# Sorted list of files in directory
simulated_list = os.listdir(options.simulated)
simulated_list = sorted(simulated_list)

derived_list = []
for fn in os.listdir(options.derived):
	derived_list.append(fn)
derived_list = sorted(derived_list)

# Iterate over sorted list.
print('Performing VEP!')
for fn in simulated_list:
	print('Working on: ' + fn) 
	
	# Determine chr number. 
	chr_num_s = fn.split('simSNPs_')[-1]
	chr_num_s = chr_num_s.replace('.vcf', '')

	# Perform commandline.
	# server with cache dir: add --dir $VEP_CACHE remove --offline
	os.system('vep --input_file '+ options.simulated + fn +' --quiet --cache --dir $VEP_CACHE --buffer 1000 --no_stats --species '+ options.species +' --format vcf --regulatory --sift b --per_gene --ccds --domains --numbers --canonical --total_length --force_overwrite --stats_file' + fn +'.html --output_file temp.vcf')
	os.system('''cat temp.vcf | awk 'BEGIN{ FS="\t"; OFS="\t"; }{ if ($1 ~ /^#/) { if ($1 ~ /^#Up/) { sub("#","",$1); print "#Chrom","Start","End",$0 } else { print } } else { split($2,a,":"); split(a[2],b,"-"); if (length(b) == 2) { print a[1],b[1],b[2],$0 } else { print a[1],b[1],b[1],$0 } }}' >> simulated_chr'''+ chr_num_s +'''_VEP-annotated.vcf''')
	

for fn in derived_list:
	print('Working on: ' + fn) 	
		
	# Determine chr number. 
	chr_num_d = fn.split('chr')[-1]
	chr_num_d = chr_num_d.replace('.vcf', '')
		
	# Perform commandline.
	os.system('vep --input_file '+ options.derived + fn +' --quiet --cache --dir $VEP_CACHE --buffer 1000 --no_stats --species '+ options.species +' --format vcf --regulatory --sift b --per_gene --ccds --domains --numbers --canonical --total_length --force_overwrite --output_file temp.vcf')
	os.system('''cat temp.vcf | awk 'BEGIN{ FS="\t"; OFS="\t"; }{ if ($1 ~ /^#/) { if ($1 ~ /^#Up/) { sub("#","",$1); print "#Chrom","Start","End",$0 } else { print } } else { split($2,a,":"); split(a[2],b,"-"); if (length(b) == 2) { print a[1],b[1],b[2],$0 } else { print a[1],b[1],b[1],$0 } }}' >> derived_chr'''+ chr_num_d +'''_VEP-annotated.vcf''')
		
# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_VEP.txt', 'x')
indication.close()
os.rename('./finished_VEP.txt', './output/finished_VEP.txt')

