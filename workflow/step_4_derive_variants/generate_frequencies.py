#!/usr/bin/env python
'''
:Author: Seyan Hu
:Date: 14-10-2022
:Usage: python <script.py> <vcf.gz file> <chr list>

Creates a frequency file per chromosome (of interest) from a vcf file.

:Example:
python scripts/generate_frequencies.py data/vcf_data/sus_scrofa.vcf.gz '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','X'
'''

# Import dependencies
import sys, os

# Assign input to variables
file_name = sys.argv[1]
chr_list = sys.argv[2].split(',')

# Loop through list of chr and perform vcftools
for chr_num in chr_list:
	print(chr_num)
	os.system('vcftools --gzvcf ' + file_name + ' --chr ' + chr_num + ' --remove-indels --non-ref-af 0.9 --max-non-ref-af 1.0 --stdout --freq > ' + chr_num + '_freq.out')


# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_generate_frequencies.txt', 'x')
indication.close()
os.rename('./finished_generate_frequencies.txt', './output/finished_generate_frequencies.txt')