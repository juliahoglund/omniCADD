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
[specied_name] -o <new path to new files>

:Example:
python sort_by_chromosome.py -p ./data/mRO/ -f mRO_mSTR_mDF_ -s sus_scrofa -o sorted


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
from collections import defaultdict

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
									<PATH TO NEW FILES> is the name of the directory that will created in which the new files will be saved")


parser.add_option_group(group)
parser.add_option("-p", "--path", dest="path", help="path to folder with processed files by mafRowOrderer", default= "./mRO")
parser.add_option("-f", "--fileprefix", dest="prefix", help="prefix of processed files (default = mRO_mStr_mDF_)", default="mRO_mStr_mDF_")
parser.add_option("-s", "--species", dest="species", help="name of species, based of which the chromosomes should be ordered", default="")
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

# Create dict for chr and corresponding blocks
maf_blocks = defaultdict(list)

# Loop through files list and add file names to correct chr
for maf_f in file_list:
	print('Processing file: {}'.format(maf_f))

	## open file if gzipped
	if maf_f.endswith(".gz"):
		gunzip(str(options.path)+str(maf_f))
		open_f = open(options.path + maf_f.replace('.gz', ''),'r')
	else: 
		open_f = open(options.path + maf_f,'r')


	ident_chr = ""
	blocks = ""
	line = open_f.readline() # header


	while line != "":
		line = open_f.readline()
		if "a " in line:
			switch = 0
			blocks = ""
			blocks += line
		elif "s " in line:
			blocks += line
			if str(options.species) in line:
				ident_chr = line.split()[1].split('.')[1]
		else:
			switch = 1
		if switch == 1:
			blocks += '\n'
			maf_blocks[ident_chr].append(blocks)
	
	os.system('gzip ' + str(options.path)+str(maf_f.replace('.gz', '')))


# Loop through keys in dict:
for chr_num, blocks in maf_blocks.items():
	if len(chrom_num) > 5:
		continue
	else:
		new_file = open("./"+options.sorted+"/"+options.species+"_chr%s.reordered.maf" %chr_num, "a")
		new_file.write("##maf version=1\n\n")
		
		for line in blocks:
			new_file.write(line)

new_file.close()

if options.clean=='yes':
	os.system('rm -r ' + str(options.path))

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_sort_by_chromosome.txt', 'x')
indication.close()
os.rename('./finished_sort_by_chromosome.txt', './output/finished_sort_by_chromosome.txt')
