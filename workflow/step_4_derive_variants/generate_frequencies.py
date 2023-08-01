#!/usr/bin/env python
'''
:Author: Seyan Hu
:Date: 14-10-2022
:Usage: python <script.py> <vcf.gz file> <chr list>

Creates a frequency file per chromosome (of interest) from a vcf file.

:Example:
python scripts/generate_frequencies.py -v data/sus_scrofa.vcf.gz -c '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X'
'''

# Import dependencies
import sys, os

# OptionParser for the directories of the input 
parser = OptionParser()
parser.add_option("-v", "--vcf", dest="vcf", help="Path to individual level vcf file", default = "reference.vcf")
parser.add_option("-c", "--chromosomes", dest="chromosomes", help="List of chromosomes to be considered when extracting freq data", default = '1,2,3,4,5')

(options, args) = parser.parse_args()


# Creates a list for the chromosomes
chr_list = options.chromosomes.split(',')


# Loop through list of chr and perform vcftools
for chr_num in chr_list:
	os.system('vcftools --gzvcf ' + options.vcf + ' --chr ' + chr_num + ' --remove-indels --non-ref-af 0.9 --max-non-ref-af 1.0 --stdout --freq > ' + chr_num + '_freq.out')

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_generate_frequencies.txt', 'x')
indication.close()

os.rename('./finished_generate_frequencies.txt', './output/finished_generate_frequencies.txt')