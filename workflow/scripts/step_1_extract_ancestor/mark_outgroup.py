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

# OptionParser for the input directory (process duplicates output), to be marked ancestor seq, and file identifier.
parser = ArgumentParser(description = __doc__)
parser.add_argument('-i', '--input',
	help = 'Alignment file to mark ancestral sequence in, can be gzip compressed', 
	required = True)
parser.add_argument('-o', '--output',
	help = 'Name of output file, will be compressed if ending in .gz',
	required = True)
parser.add_argument('-a', '--ancestor',
	help = 'The name used to mark the outgroup sequence as ancestor',
	required = True)
parser.add_argument('--sp1-label',
	help = 'alignment label/name of the species of interest e.g. sus_scrofa', 
	required = True)

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
	# reformat input to fit maf format
	name_sp1 = 's ' + args.sp1_label + '.'
	name_ancestor = 's ' + args.ancestor + '.'

	try:
		alignment_file = gzip.open(args.input, "rt") \
			if args.input.endswith('.gz') else open(args.input, "r")
	except IOError as e:
		print(f"Error opening input file: {e}")
		return

	try:
		outfile = gzip.open(args.output, "wt", compresslevel=1) \
			if args.output.endswith('.gz') else open(args.output, "w")
	except IOError as e:
		print(f"Error opening output file: {e}")
		alignment_file.close()
		return

	for lines in alignment_file:
		lines = lines.strip()
		if not lines:
			continue
		if lines.startswith(name_sp1):
			# Write out given species
			outfile.write('\t'.join(lines.split()) + '\n')
		elif lines.startswith(name_ancestor):
			newline = lines.strip().split('.')[1]
			outfile.write(f's\tAncestor_{args.ancestor}.' + '\t'.join(newline.split()) + '\n')
		elif lines.startswith('##maf'):
			continue
		else:
			outfile.write('\t'.join(lines.split()) + '\n')

	outfile.close()
	alignment_file.close()

if __name__ == '__main__':
	main(parser.parse_args())

