#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
:Original Author: Christian Gross
:Contact: cgross@tudelft.nl
:Date: 18-11-2017

If the user does not want to extract a reconstructed ancestor, or if there is non
in the given MSA file, this script will extract and mark a chosen species instead.
By adding "Ancestor_" before the name of the species in the msa, 
The rest of the pipeline will recognize it as the ancestral sequence in the rest of the pipeline

:Updating and reworking
:Author: Seyan Hu
:Date: 26-9-2022

:Extension and modification
:Author: Julia Höglund
:Contact: julia.hoglund@su.se
:Date: 27-2-2023

:Usage: python mark_outgroup.py -p <path to files> -a [ancestor (outgroup to be marked)] -i [file identifier] -f [scientific name]

:Example:
python mark_outgroup.py -p ./ -a bos_taurus -i 43_mammals.epo -f sus_scrofa
"""
# Import dependencies
import sys, os, argparse
from optparse import OptionParser, OptionGroup
import gzip
from Bio import Phylo
from Bio import AlignIO
from io import StringIO
from sh import gunzip

# OptionParser for the input directory (process duplicates output), to be marked ancestor seq, and file identifier.
parser = OptionParser()
parser.add_option("-p", "--path", dest="path", help="path to folder with alignment files", default= "./")
parser.add_option("-a", "--ancestor", dest="ancestor", help="sequence which should be marked as ancestror", default= "")
parser.add_option("-i", "--identifier", dest="identifier", help="identifier to get the correct files; (eg. chr1[.maf])", default="")
parser.add_option("-f", "--full", dest="name", help="full scientific name of species of interest, separated with an underscore; (eg. homo_sapiens)", default="")

if len(sys.argv)==1:
   	parser.print_help()
   	parser.exit()
(options, args) = parser.parse_args()

# CHECK input
if not os.path.exists(options.path):
 	sys.stderr.write("Invalid path to given msa directory.\n")
 	sys.exit()
# Making sure that the path ends with '/'
if (not options.path.endswith('/')):
	options.path = options.path+'/'

if options.identifier == '':
 	sys.stderr.write("Found no files with given identifier (prefix).\n")
 	sys.exit()
if options.name == '':
  	sys.stderr.write("No full species name given!\n")
  	sys.exit()

# reformat input to fit maf format
name_sp1 = 's ' + options.name + '.'
name_ancestor = 's ' + options.ancestor + '.'


# Creates a list of all files.
file_list = os.listdir(options.path)
# creates a list of files to be analysed
processed_file_list = []
for file in file_list:
	if (options.identifier in file) and (not 'README' in file):
		processed_file_list.append(file)

## if path not exist makedir marked else not do it
if not os.path.isdir('marked'):
    os.system('mkdir marked')

for file in processed_file_list:
	#####################################################
	# Analysing trees in maf
	#####################################################

	print('Marking ancestor in file: {}'.format(file))

	outfile = open('marked/marked_'+file.replace('.gz', ''), "w")

	if file.endswith(".gz"):
		gunzip(options.path+file)
	alignment_file = open(options.path+file.replace('.gz', ''), "r")
	for lines in alignment_file:
		if lines.startswith(name_sp1):
			# Write out given species
			outfile.write(lines)
		elif lines.startswith(name_ancestor):
			newline = lines.strip().split('.')[1]
			outfile.write('s Ancestor_' + options.ancestor + '.' + newline)
		else:
			outfile.write(lines)
	os.system('gzip marked/marked_' + file.replace('.gz', ''))
	outfile.close()
	os.system('gzip ' + str(options.path+file.replace('.gz', '')))
	alignment_file.close()

## if path not exist makedir output else not do it
if not os.path.isdir('output'):
    os.system('mkdir output')

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_mark_outgroup.txt', 'x')
indication.close()
os.rename('./finished_mark_outgroup.txt', './output/finished_mark_outgroup.txt')
