#!/usr/bin/env python
'''
:Author: Seyan Hu
:Date: 2-11-2022
:Extenion and modification: Julia Höglund
:Date: 02-05-2023

:Usage: python <script.py> -a <ancestor seq folder> -r <reference seq folder> -c <chromosome number> -p <ancestor seq file prefix> -r <ref seq file prefix> -s <ref species name>

:Example:
python scripts/wrapper_create_parameters.py -a ./output/extracted_ancestor -r ./genome/ -c 1,2,3,4,5 -p Ancestor_Human_ -r Homo_sapiens_ref_ -s homo_sapiens

This script is a wrapper for create_parameters.py. 
Since that script can only take one input file, 
this script is used as a wrapper to loop through multiple input files.

'''

# Import dependencies
import sys, os
from optparse import OptionParser


# OptionParser for the directories of the input 
parser = OptionParser()
parser.add_option("-a", "--ancestor-path", dest="an_p", help="Path to ancestor folder", default= './output/extracted_ancestor/')
parser.add_option("-g", "--genome-path", dest="r_p", help="Path to (reference) genome folder", default= './genome/')
parser.add_option("-c", "--chromosomes", dest="chr", help="List of chromosomes (that have an ancestral sequence)", default= "1")
parser.add_option("-p", "--prefix", dest="prefix", help="prefix of ancestral files", default= "Ancestor_")
parser.add_option("-s", "--species", dest="species", help="name of species of interest",default= 'homo_sapiens')
parser.add_option("-r", "--reference", dest="reference", help="prefix of reference species files", default= "Sus_scrofa_ref_")


(options, args) = parser.parse_args()


# Creates a list for the chromosomes
chr_list = options.chr.split(',')

# Loops through each chromosome and performs the commandline
for chr_number in chr_list:
	os.system('python scripts/create_parameters.py -a ' + options.an_p + options.prefix + chr_number + '_chr' + chr_number + '.fa' + ' -r ' + options.r_p + options.reference + chr_number + '.fa' + ' -c ' + chr_number + ' -o ' + options.species + '_chr' + chr_number + '.log')

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_substitution_calc.txt', 'x')
indication.close
os.rename('./finished_substitution_calc.txt', './output/finished_substitution_calc.txt')


