#!/usr/bin/env python
'''
:Author: Seyan Hu
:Date: 29-9-2022
:Extension and modification: Julia Höglund
:Date: 22-3-2023

:Usage: python wrapper_extract_ancestor.py -p <path to input files> -a [ancestor name] -s [scientific name of ref. species] -f [prefix of maf files] -g <path to gen ancestor script>
:Example:
python wrapper_extract_ancestor.py -p ./output -a Ancestor_ancestral_sequence -s homo_sapiens -f rmOS_ -g ./scripts -i genome.fai

This script is a wrapper for extract_ancestor.py. 
Since that script can only take one input file, 
this script is used as a wrapper to loop through multiple input files.

It is used to run this command line, looping through multiple files: 
python scripts/extract_ancestor.py output/fileN.maf Ancestor_ancestral_seq homo_sapiens 19
'''


# Import dependencies
import sys, os
from optparse import OptionParser

# OptionParser for the directories of the input 
parser = OptionParser()
parser.add_option("-p", "--path", dest="path", help="Path to maf alignment file.", default='./forwardStrandOnly') 
parser.add_option("-a", "--ancestor", dest="ancestor", help="Ancestor identifier of the ancestral sequence which should be retrieved.", default= 'Ancestor_ancestral_sequence')
parser.add_option("-s", "--species", dest="species", help="Scientific name of the species of interest.", default='')
parser.add_option("-f", "--fileprefix", dest="prefix", help="Prefix of the maf files from which the ancestor will be extracted", default='rmOS_rmSP_mS_')
parser.add_option("-g", "--generate", dest="generate", help="path to where the 'generate_ancestor.py' script is located", default='./scripts/')
parser.add_option("-i", "--fa-index", dest="index", help="full path to the indexfile to the reference genome of the species of interest", default='./genome.fai')

(options, args) = parser.parse_args()

# CHECK input
if not os.path.isdir(options.path):
	sys.stderr.write("Found no directory with that name!\t{}\n".format(options.path))
	sys.exit()
if not os.path.isdir(options.generate):
	sys.stderr.write("Found no directory with that name!\t{}\n".format(options.generate))
	sys.exit()
if options.species == '':
	sys.stderr.write("No species name given.\n")
	sys.exit()

# Checking if the path ends with '/' 
if (not options.path.endswith('/')):
	options.path = options.path+'/'

# Makes a list from the files in the directory
file_list = []
for maf_file in os.listdir(options.path):
	if maf_file.startswith(options.prefix):
		file_list.append(maf_file)

if bool(file_list) == False:
	sys.stderr.write("Found no files with the given prefix, {}\n".format(options.prefix))
	sys.exit()

# Loops through each file in given path and performs the commandline
for inp in file_list:
	os.system('python3 ' + options.generate + '/extract_ancestor.py ' + options.path + inp + ' ' + options.ancestor + ' ' + options.species + ' ' + options.index)
	file = options.path + '_'.join(inp.split('_')[-3:])
	os.system('mv ' + str(options.path) + str(inp) + ' ' + str(file))

os.system('mv ' + str(options.path) + ' processedMAFfiles')
print("Last directory in pipeline have been renamed from {} to processedMAFfiles.".format(options.path))
print("Fully preprocessed files have been renamed to {}_chrN.maf.".format(options.species))
print("Step 1; Extract ancestral sequence done.\n")

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_extract_ancestor.txt', 'x')
indication.close
os.rename('./finished_extract_ancestor.txt', './output/finished_extract_ancestor.txt')

