#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
:Author: Seyan Hu
:Date: 29-9-2022
:Extension and modification: Julia Höglund
:Date: 2023-03-17

This script takes as input the path to the folder containing  maf files.
Each file in each chromosome will be read and then sorted to the correct chromosome.
This is performed to create one maf file per chromosome,
since there are multiple maf files per chromosome currently.

:Usage: python sort_by_chromosome.py -p <path to processed files> -f <prefix of the processed files> -s 
[specied_name] -k [chr1,chr2,chrn...] -o <new path to new files>

:Example:
python sort_by_chromosome.py -p ./data/mRO/ -f mRO_mSTR_mDF_ -s sus_scrofa -k 
'1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X,Y,Other' -o sorted


"""

# Import dependencies
import sys, os
from optparse import OptionParser, OptionGroup
import gzip
from sh import gunzip
from Bio import Phylo
from Bio import AlignIO
from io import StringIO
from Bio.Align import MultipleSeqAlignment

# OptionParser for the input directory (process duplicates output), to be marked ancestor seq, and file identifier.
parser = OptionParser()
group = OptionGroup(parser, "Usage (extended)", "This script takes as input the maf files where the \
            ancestor of choice has been marked, the reference has been extracted, \
            and the file have been processed by mafTools. \
				The input is given like this: \
						python sort_by_chromosome.py -p <path to processed files> -f <prefix of the processed files> -s [species_name] -k [chr1,chr2,chrn...] -o <path to new files>\
Where: \
								<PATH TO PROCESSED FILES> is the path to the files that has been processed in the previous step by mafTools, 'apply_mafTools' \
		<PREFIX OF THE PROCESSED FILES> is the prefix of the processed files \
[SPECIES NAME] iis the name of the species, based on which the chromosomes and alignments should be ordered \
			[CHR1,CHR2,CHRn...] is the chromosomes to be considered, specified one by one \
									<PATH TO NEW FILES> is the name of the directory that will created in which the new files will be saved")


### LÄGG TILL ATT GÖRA EN DIR FÖR OUTPUT ATT SPARA ALLT I SOM MATCHAR NÄSTA
parser.add_option_group(group)
parser.add_option("-p", "--path", dest="path", help="path to folder with processed files by mafRowOrderer", default= "./mRO")
parser.add_option("-f", "--fileprefix", dest="prefix", help="prefix of processed files (default = mRO_mStr_mDF_)", default="mRO_mStr_mDF_")
parser.add_option("-s", "--species", dest="species", help="name of species, based of which the chromosomes should be ordered", default="")
parser.add_option("-k", "--chromosomes", dest="chromosomes", help="chromosomes that should be considered in the ordering, including possible sex chromosomes and 'Other'", default="1,X,Y")
parser.add_option("-o", "--ordered", dest="sorted", help="directory within which the new output will be saved", default="sorted")
# clean intermediate files
parser.add_option("-c", "--clean", dest="clean", help="remove intermediate files (from previous step) after computation [yes/no]", default="no")

if len(sys.argv)==1:
    parser.print_help()
    parser.exit()
(options, args) = parser.parse_args()

# CHECK input
if not os.path.isdir(options.path):
    sys.stderr.write("Found no directory with that name, {} \n".format(options.path))
    sys.exit()
if options.species == '':
    sys.stderr.write("No species name given.\n")
    sys.exit()
if options.chromosomes == '':
    sys.stderr.write("Forgot to specify chromosomes.\n")
    sys.exit()

# Checking if the path ends with '/'
if (not options.path.endswith('/')):
	options.path = options.path+'/'

## if path not exist makedir sorted else not do it
if not os.path.isdir(options.sorted):
    os.system('mkdir ' + str(options.sorted))

# Makes a list from the files in the directory and remove 'other' from files
file_list = []
for file in os.listdir(options.path):
	if file.startswith(options.prefix):
		file_list.append(file)

if bool(file_list) == False:
	sys.stderr.write("Found no files with the given prefix.\n")
	sys.exit()

# Create dict for chr and corresponding files
output_files_chr_dict = {}
chromosomes = options.chromosomes.split(',')
for item in chromosomes:
	output_files_chr_dict[item] = []

# Loop through files list and add file names to correct chr
for maf_f in file_list:
	print('Processing file: {}'.format(maf_f))

	## open file if gzipped
	if maf_f.endswith(".gz"):
		gunzip(str(options.path)+str(maf_f))
		open_f = open(options.path + maf_f.replace('.gz', ''),'r')
	else: 
		open_f = open(options.path + maf_f,'r')

	content_f = open_f.readlines()

	for lines in content_f:
		if str(options.species) in lines:
			if "AEMK" in lines:
				break
			correct_l = lines
			ident_chr = correct_l.split()[1].split('.')[1]
			output_files_chr_dict[ident_chr].append(maf_f)
			break

# Loop through keys in dict:
for chr_num, file_names in output_files_chr_dict.items():
	new_file = open("./"+options.sorted+"/"+options.species+"_chr%s.reordered.maf" %chr_num, "a")
	new_file.write("##maf version=1\n\n")
	for fn in file_names:
		open_inp = open(options.path + fn,'r')

		for line in open_inp.readlines():
			new_file.write(line)

		os.system('gzip ' + str(options.path)+str(fn))

new_file.close()

if options.clean=='yes':
	os.system('rm -r ' + str(options.path))

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_sort_by_chromosome.txt', 'x')
indication.close
os.rename('./finished_sort_by_chromosome.txt', './output/finished_sort_by_chromosome.txt')

