#!/usr/bin/env python
'''
:Author: Seyan Hu
:Date: 11-10-2022
:Extension and modification: Julia Höglund
:Date: 20-3-2023
:Usage: python remove_species.py -p <path to input files> -f <file prefix> -s <species> -r <new path to new files>

This script removes unwanted species from maf files.
This script should be performed after mafSorter in the pipeline. 

:Example:
python scripts/remove_species.py -p output/ -f mS_ -s mus_musculus pruned
'''

# also make sure the name is like understandable

# Import dependencies
import sys, os
from optparse import OptionParser
import gzip
from sh import gunzip

# Option parser for path
parser = OptionParser()
parser.add_option("-p", "--path", dest="path", help="path to folder with mafSorted_EPO files.", default= './sortedMSA') 
parser.add_option("-f", "--fileprefix", dest="prefix", help="prefix of the processed files (default mS_)", default='mS_')
parser.add_option("-s", "--species", dest="species", help="species of interest (i.e. reference species)", default='')
parser.add_option("-r", "--removed", dest="removed", help="directory within which the new output will be saved", default='pruned')
parser.add_option("-c", "--clean", dest="clean", help="remove intermediate files after computation [yes/no]", default="no")

if len(sys.argv)==1:
    parser.print_help()
    parser.exit()
(options, args) = parser.parse_args()

# CHECK input
if not os.path.isdir(options.path):
    sys.stderr.write("Found no directory with that name, {}\n".format(options.path))
    sys.exit()
if options.species == '':
    sys.stderr.write("No species name given.\n")
    sys.exit()

# Checking if the path ends with '/'
if (not options.path.endswith('/')):
	options.path = options.path+'/'

## if path not exist makedir pruned else not do it
if not os.path.isdir(options.removed):
    os.system('mkdir ' + str(options.removed))

# Makes a list from the files in the directory
file_list = []
for maf_file in os.listdir(options.path):
	if maf_file.startswith(options.prefix):
		file_list.append(maf_file)

if bool(file_list) == False:
	sys.stderr.write("Found no files with the given prefix, {}\n".format(options.prefix))
	sys.exit()

# Loop through files in file_list
for file_name in file_list:
	open_f = open(options.path + file_name,'r')
	content_f = open_f.readlines()
	
	print("Processing file: {}".format(file_name))

	list_lines = []
	# Loop through lines in file, 
	# add lines to list that doesn't contain species with similar name
	for line in content_f:
		if not line.startswith('s ' + options.species + '_'):
			list_lines.append(line)	

	# Open new file
	# Write lines in list to new file.
	new_f = open(options.removed + '/rmSP_' + file_name, 'w')
	new_f.writelines(list_lines)
	new_f.close

if options.clean=='yes':
	os.system('rm -r ' + str(options.path))

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_removing_unwanted_species.txt', 'x')
indication.close
os.rename('./finished_removing_unwanted_species.txt', './output/finished_removing_unwanted_species.txt')


