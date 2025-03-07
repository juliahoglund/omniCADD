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

Edited by Job van Schipstal, 15-9-2023. Using argparse for input
:Extension and modification
instead. Cleaned up warnings for potentially undefined objects. 
Formatted file to satisfy pylint.

:Extension and modification
:Author: Julia Höglund
:Contact: julia.hoglund@su.se
:Date: 20-10-2023

:Usage: python mark_ancestor.py -p <path to files> -a [ancestor] -i [file
    identifier] -s [species1,species2] -c [abbreviation1,abbreviation2] -f
    [scientific name]

:Example:
python mark_ancestor.py -p ./ -a Pig_Cow_Sheep -i 43_mammals.epo -s Pig,Cow -c Sscr,Btau -f sus_scrofa
"""

# Import dependencies
from argparse import ArgumentParser
import gzip
from Bio import Phylo
from io import StringIO
import logging
import sys
import os

# Set up logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

parser = ArgumentParser(description=__doc__)
parser.add_argument('-i', '--input',
	help='Alignment file to mark ancestral sequence in, can be gzip compressed', 
	required=True)
parser.add_argument('-o', '--output',
	help='Name of output file, will be compressed if ending with .gz',
	required=True)
parser.add_argument('-a', '--ancestor',
	help='The name used to mark the ancestral sequence, eg. Pig_Cow',
	required=True)
parser.add_argument('-l', '--logfile',
	help='Log file (default: mark_ancestor.log)',
	default="mark_ancestor.log")
parser.add_argument('--sp1-label',
	help='alignment label/name of the species of interest e.g. sus_scrofa', 
	required=True)
parser.add_argument('--sp1-ab',
	help='abbrevation of species of interest as in the tree e.g. Sscr', 
	required=True)
parser.add_argument('--sp2-ab',
	help='abbrevation of species 2 as in the tree e.g. Btau',
	required=True)

def hierarchy(found_clades, args):
	"""
	This function takes in a list of Bio.Phylo.BaseTree.Clade objects,
	which represent the terminal taxa of the investigated tree. The function
	iterates over these taxa to identify if the given species is available and
	returns the taxa which should be used for the identification of the last
	common ancestor of the given species.
	"""
	# Abbreviation of the scientific name of species 1 and 2 with '_' at end
	ab_sp1 = args.sp1_ab + '_'
	ab_sp2 = args.sp2_ab + '_'
	sp1 = None
	sp2 = None

	for f_clade in found_clades:
		if (ab_sp1 in f_clade.name) and (not ab_sp1 + '_AEMK' in f_clade.name):
			sp1 = f_clade.name
		elif ab_sp2 in f_clade.name:
			sp2 = f_clade.name

	if sp1 and sp2:
		return [sp1, sp2]
	return []

def main(args):
	"""
	Main function. Opens relevant file. Iterates through lines in maf file.
	Iterate until I hit '# tree:' and read in the newick string and identify
	if there is a common ancestor. If yes, write out the ancestral
	identifier and the information to which chromosome it was aligned,
	to save the ancestral sequence according to the chromosome it belongs,
	this will ease the aligning effort, since it will only have to align to
	the chromosome it belongs to and will not have to take all chromosomes
	into account.
	:param args: argparse parameters
	:return: None
	"""

	try:
		with gzip.open(args.input, "rt") if args.input.endswith('.gz') else open(args.input, "r") as alignment_file:
			with gzip.open(args.output, "wt", compresslevel=1) if args.output.endswith('.gz') else open(args.output, "w") as outfile:
				with open(args.logfile, 'w') as log_s:
					log_s.write(f"Logfile from mark_ancestor.py\nScript arguments:\n{args}\n")

					ancestor_id = ''
					i = 0
					elem = None
					ancestor_seq = None

					for lines in alignment_file:
						to_write = True
						# Checking if tree contains information regarding ancestor
						if '# tree:' in lines:
							ancestor_seq = None
							elem = None
							i += 1
							lines = lines.strip()
							to_write = False

							# Processing the tree to identify the ID of the last common ancestor
							newick = lines.strip().split(' ')[-1]
							tree = Phylo.read(StringIO(newick), "newick")

							clades = hierarchy(tree.get_terminals(), args)

							if len(clades) == 0:
								ancestor_id = 'not_available'
								log_s.write(str(i) + '\t' + ancestor_id + '\t' + 'NA' + '\n')
							else:
								ancestor_id = '_'.join(tree.common_ancestor(clades).name.split('_')[1:4])
								log_s.write(str(i) + '\t' + ancestor_id + '\t' + args.ancestor + '\n')
						elif lines.startswith(f's {args.sp1_label}.'):
							# Stores coords for ancestor
							elem = lines.strip().split()[0:-1]
						elif lines.startswith('s ancestral_sequences.' + ancestor_id):
							to_write = False
							ancestor_seq = lines.strip().split()[-1]

						# Write renamed/annotated ancestral alignment if fields found.
						if ancestor_seq and elem:
							outfile.write('s ' + args.ancestor + '.' +	
								elem[1].split('.')[1] + '\t' + '\t'.join(elem[2:]) + '\t' + ancestor_seq + '\n')
							ancestor_seq, elem = None, None

						# Write original line if not tree or ancestral sequence
						if to_write:
							outfile.write(lines)
	except Exception as e:
		logging.error(f"Error processing files: {e}")
		sys.exit(1)

if __name__ == '__main__':
	main(parser.parse_args())
