#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
:Original Author: Christian Gross
:Contact: cgross@tudelft.nl
:Date: 18-11-2017

Identifying the last common ancestor and marking them.
This is necessary, thus they are not removed when applied to mafTools to remove duplicates.
Furthermore, this script is setting the coordinates equal to the wanted species,
so it will be flipped if your specific species is also flipped.

:Updating and reworking
:Author: Seyan Hu
:Date: 26-9-2022

:Extension and modification
:Author: Julia Höglund
:Contact: julia.hoglund@su.se
:Date: 27-2-2023

:Usage: python mark_ancestor.py -p <path to files> -a [ancestor] -i [file
    identifier] -s [species1,species2] -c [abbreviation1,abbreviation2] -f
    [scientific name]

:Example:
python mark_ancestor.py -p ./ -a Pig_Cow_Sheep -i 43_mammals.epo -s Pig,Cow -c Sscr,Btau -f sus_scrofa
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
group = OptionGroup(parser, "Usage (extended)", "This script takes as input the desired species from which \
the most recent common ancestor should be found and marked. \
			The input is given like this: \
						python mark_ancestor.py -p <path to files> -a [ancestor] -i [file identifier] -s [species1,species2] -c [abbreviation1,abbreviation2] -f [scientific name] \
							Where: \
									<PATH TO FILES> is the path to the msa files (maf, gzipped or standard) to be marked \
									<ANCESTOR> is the name of the ancestor that should be marked \
									<FILE IDENTIFIER> is the identifier (prefix) of the maf files \
<SPECIES1,SPECIES2> is the names of the species (separated by a comma, no space) from which the common ancestor will be marked \
<ABBREVIATION1,ABBREVIATION2> is the abbreviated species names, as stated in the .maf file \
						<SCIENTIFIC NAME> is the full scientific name of the species of interest \
								Based on this, a log file together with (new) marked .maf files will be created.")

parser.add_option_group(group)
parser.add_option("-p", "--path", dest="path", help="path to folder with alignment files", default= "./")
parser.add_option("-a", "--ancestor", dest="ancestor", help="ancestor sequence which is should be marked", default= "Ancestral_Sequence")
parser.add_option("-i", "--identifier", dest="identifier", help="identifier to get the correct files; (eg. chr1[.maf])", default="")
parser.add_option("-s", "--species", dest="species", help="name of species one (1) and species two (2); (eg. Human,Chimp)", default="")
parser.add_option("-c", "--compressed", dest="abbreviation", help="abbreviation of the given species as stated in the .maf file; (eg. Hsap,Ptro)", default="")
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
if options.species == '':
 	sys.stderr.write("No species names given!\n")
 	sys.exit()
if options.abbreviation == '':
  	sys.stderr.write("No abbreviated species names given!\n")
  	sys.exit()
if options.name == '':
  	sys.stderr.write("No full species name given!\n")
  	sys.exit()

# reformat input to fit maf format
ab_sp1 = options.abbreviation.split(',')[0] + '_'
ab_sp2 = options.abbreviation.split(',')[1] + '_'
sp1 = options.species.split(',')[0]
sp2 = options.species.split(',')[1]
name_sp1 = 's ' + options.name + '.'

# This function takes in a list of Bio.Phylo.BaseTree.Clade objects,
# which represent the terminal taxa of the investigated tree.
# The function iterates over these taxa to identify if the given species is available
# and returns the taxa which should be used
# for the identification of the last common ancestor of the given species
# and the name which should be given to that ancestor.
def hierarchy(found_clades,required_ancestor):
	binary_dict = {}
	binary_dict[sp1] = ''
	binary_dict[sp2] = ''


	for f_clade in found_clades:
		if (ab_sp1 in f_clade.name) and (not ab_sp1 + '_AEMK' in f_clade.name):
			binary_dict[sp1] = f_clade.name
		elif ab_sp2 in f_clade.name:
			binary_dict[sp2] = f_clade.name

	# By applying this order of return functions it introduces a hierarchy
	if (not binary_dict[sp1]=='') and (not binary_dict[sp2]=='') and (required_ancestor==options.ancestor):
		return [binary_dict[sp1],binary_dict[sp2]],'Ancestor_' + options.ancestor
	else:
		return [],None

# Creates a list of all files.
file_list = os.listdir(options.path)
# creates a list of files to be analysed
processed_file_list = []
for file in file_list:
	if (options.identifier in file) and (not 'README' in file):
		processed_file_list.append(file)
log_file = open(options.path+'%s.log'%(options.ancestor),'w')

## if path not exist makedir marked else not do it
if not os.path.isdir('marked'):
    os.system('mkdir marked')

for file in processed_file_list:
	#####################################################
	# Analysing trees in maf
	#####################################################


	# Iterate over the list until it hits '# tree:'
	# and read in the newick string and identify if there is a common anecestor
	# If yes, write out the ancestral identifier and
	# the information to which chromosome it was aligned,
	# to save the ancestral sequence according to the chromosome it belongs,
	# this will ease the aligning effort,
	# since it only have to align it to the chromosome it belongs to
	# and do not have to take all chromosomes into account.
	print('Marking ancestor in file: {}'.format(file))

	outfile = open('marked/marked_'+file.replace('.gz', ''), "w")
	log_file.write('./'+file+'\n')
	ancestor_id = ''
	i = 0
	if file.endswith(".gz"):
		gunzip(options.path+file)
	alignment_file = open(options.path+file.replace('.gz', ''), "r")
	for lines in alignment_file:
		# Checking if tree contains information regarding ancestor
		if '# tree:' in lines:
			i += 1
			# Tree is always comes before alignment block,
			# therefore, the trees can be used to reset variables which may
			# be filled by some values below
			switch = False
			ancestor_seq = ''
			elem = []
			# Processing the tree to idenfy the ID of the last common ancestor
			newick = lines.strip().split(' ')[-1]
			handles = StringIO(newick)
			tree = Phylo.read(handles, "newick")
			clades,new_ancestor_name = hierarchy(tree.get_terminals(),options.ancestor)
			if len(clades)==0:
				ancestor_id = 'not_available'
				log_file.write(str(i)+'\t'+ancestor_id+'\t'+'NA'+'\n')
			else:
				ancestor_id = '_'.join(tree.common_ancestor(clades).name.split('_')[1:4])
				log_file.write(str(i)+'\t'+ancestor_id+'\t'+new_ancestor_name+'\n')
				# The following two code blocks are changing
				# depending on if the ancestor sequence appears before the sequence
				# or the other way around.
		if lines.startswith(name_sp1):
			# Stores coords for ancestor
			elem = lines.strip().split()[0:-1]
			# Write out given species
			outfile.write(lines)
			if ancestor_seq == '':
				switch = True
				continue
			else:
				# Do some changes and write out ancestor
				outfile.write('s '+new_ancestor_name+'.'+elem[1].split('.')[1]+'\t'+'\t'.join(elem[2:])+'\t'+ancestor_seq+'\n')
				switch = False
				continue
		# This if statement is never entered if there is no name
		# of given species in the alignment block
		if lines.startswith('s ancestral_sequences.'+ancestor_id):
			if switch==True:
				# If coords are already known, do some changes to ancestral and write out ancestral
				outfile.write('s '+new_ancestor_name+'.'+elem[1].split('.')[1]+'\t'+'\t'.join(elem[2:])+'\t'+lines.strip().split()[-1]+'\n')
				switch = False
				continue
			else:
				# Store ancestor if given species has not been called up yet
				ancestor_seq = lines.strip().split()[-1]
				continue
		outfile.write(lines)
	log_file.flush()
	os.system('gzip marked/marked_' + file.replace('.gz', ''))
	outfile.close()
	os.system('gzip ' + str(options.path+file.replace('.gz', '')))
	alignment_file.close()

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_mark_ancestor.txt', 'x')
indication.close()
os.rename('./finished_mark_ancestor.txt', './output/finished_mark_ancestor.txt')
