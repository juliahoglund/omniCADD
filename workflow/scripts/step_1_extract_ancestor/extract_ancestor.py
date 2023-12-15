#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script takes as input the preprocessed maf file, the ancestor sequence
of interest and returns per .maf file (per chromosome) the ancestor sequence
in fasta format.

:Author: Christian Gross
:Contact: cgross@tudelft.nl
:Date: 01.02.2018

:Edited by: Seyan Hu :Date: 29-9-2022
Added: 1) It removes ancestral sequences that are not of the same length as
stated in their annotation. 2) It checks whenever the start position of the
to be added sub sequence is located in the previous added sub sequence. If
true, this sequence will be skipped.

:Edited by: Seyan Hu
:Date: 29-9-2022
:Modified by: Julia Höglund and Julia Beets
:Date 22-3-2023

Added: 
1) It removes ancestral sequences that are not of the same length as stated in their annotation. 
2) It checks whenever the start position of the to be added sub sequence is located in the previous added sub sequence.
	If true, this sequence will be skipped. 

:Edited by: Job van Schipstal :Date: 19-9-2023
:Example Usage: python scripts/gen_ancestor_seq.py
sorted_chr19.maf ancestor_chr19.fa Ancestor_Mouse_Rat mus_musculus

Reformatted script according to PEP guidelines, documented functions, extracted
get_alignment_generator and get_ancestral_sequence functions from main code.
"""

import os
import sys
import gzip
from argparse import ArgumentParser
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
	help = "maf alignment file from which the ancestral sequence is to be extracted. It is expected that it has been filtered and sorted.", 
	required = True)
parser.add_argument("-o", "--output", 
	help = "Ancestral sequence output file, fasta format", 
	required = True)
parser.add_argument("-a", "--ancestor", 
	help = "Name/label of the ancestral sequence", 
	required = True)
parser.add_argument("-n", "--interest", 
	help="Name/label of the species of interest",
	required = True)
parser.add_argument("-r", "--reference", 
	help="Name of the genomewide reference species fasta sequence",
	required = True)

args = parser.parse_args()

def ancestor_finder(alignment_instance, anc_identifer):
	"""
	A function that checks if the ancestor sequence of interest
	is available in alignment block and returns seq, otherwise empty string
	:param alignment_instance: AlignIO alignment block to search in
	:param anc_identifer: prefix of ancestral id to look for
	:return: entry with ancestral id, None if not found
	"""
	for idx, entrance in enumerate(alignment_instance):
		if entrance.id.startswith(anc_identifer):
			return entrance

	return None


def sequence_processing(g_seq, ancestor_seq):
	"""
	Function that removes gaps that are present in both the sequence of the
	given species and its ancestor. This is done to find the real genomic
	location of the ancestor to the reference. Otherwise it would result in
	wrong positions when deriving variants.
	:param g_seq: seqRecord of reference/target sequence
	:param ancestor_seq: seqRecord of ancestor sequence
	:return: Processed ancestral sequence, type seqRecord or str 'Removed'
	"""
	# Initializing two empty SeqRecord objects which will be filled with the
	# processed given species and ancestor sequences respectively. This way
	# the true length of the ancestor can be identified and it will fill up
	# all the remaining positions with gaps.
	anc_out_seq = SeqIO.SeqRecord(id=ancestor_seq.id, seq=Seq(''),
		name=ancestor_seq.name, annotations={
		'start': ancestor_seq.annotations['start'],
		'srcSize': ancestor_seq.annotations['srcSize'],
		'strand': ancestor_seq.annotations['strand'],
		'size': ancestor_seq.annotations['size']})

	for (g_char, anc_char) in zip(g_seq, ancestor_seq):
		if g_char != '-':
			anc_out_seq.seq = anc_out_seq.seq + Seq(anc_char)

	# Checks if the length is the same as the size in annotations
	if len(anc_out_seq) != anc_out_seq.annotations['size']:
		sys.stderr.write(
			'The length of the processed ancestor sequence differs from the '
			'stated size: \n Annotations : {0} \n length : {1}\n'.format(anc_out_seq.annotations, len(anc_out_seq)))
		return 'Removed'
	return anc_out_seq


def get_alignment_generator(input_f):
	"""
	Reads input alignment using AlignIO. supports either a
	.maf file or .maf.gz if compressed
	:param input_f: path to input file, maf or maf.gz format.
	:return: AlignIO parser on the input file
	"""
	if not os.path.isfile(input_f):
		sys.exit('Program stopped prematurely\n Parsed path does not point to a '
			'file. The parsed file is following : {0}'.format(input_f))

	handle = gzip.open(input_f, "rt") if str(input_f).endswith(".gz") \
		else open(input_f, "r")
	return AlignIO.parse(handle, "maf")


# get chromosomes to use
if 'chr' in args.input:
	chr_num = args.input.split('chr')
	chromosome = chr_num[1].split('.')[0]
else:
	chr_num = 'N'

# get total chromosome length (in bases)
chr_lengths = defaultdict(dict)

with open(args.reference) as f:
	lines = f.read().splitlines()

lines = [ x for x in lines if len(x.split('\t')[0]) < 6 ]

for line in lines:
	chrom = line.split()[0]
	length = line.split()[1]
	chr_lengths[chrom] = length

# Init dictionary with sequence start location as key and the sequence as values.
# After iterating over the alignments it will go through this dictionary 
# and create the final ancestor sequence by filling up missing positions with gaps.
anc_pos_seq_dict = {}

print("## Reading MAF file")
# Loop through the file and create a dict for the ancestral seq. 
for alignment in get_alignment_generator(args.input):
	sp_seq = alignment[0]	
	
	# Search for the ancestor sequence of interest in alignment block and checks if the seq is of same size as in the annotation
	anc_seq = ancestor_finder(alignment,args.ancestor)
	ref_seq = ancestor_finder(alignment,args.interest)
	
	if (not isinstance(anc_seq,SeqIO.SeqRecord)) or (str(args.interest) +'.' not in sp_seq.id):
		if (isinstance(anc_seq,SeqIO.SeqRecord)) and (str(args.interest) +'.' not in sp_seq.id):
			print('Something went wrong in the processing. The alignment block contains an ancestor sequence but not '+ str(path_to_input) +' at first position. Path and sequence details %s \t %s' %(path_to_input,anc_seq.annotations))
			print('This alignment block will be skipped.\n')
			continue
		else:
			continue
	else:
		#ADDED TO MAKE POS REF GENOME, AND SAME START COORDINATES
		anc_seq_processed = sequence_processing(sp_seq, anc_seq)
		ref_start = ref_seq.annotations['start']
		if ref_start in anc_pos_seq_dict:
			if len(anc_seq_processed) > len(anc_pos_seq_dict[ref_start]):
				anc_pos_seq_dict[ref_start] = anc_seq_processed
			else:
				continue
		else:
			anc_pos_seq_dict[ref_start] = anc_seq_processed

# Read alignment using AlignIO, extract ancestral sequences and put them in
# dict start_pos=sequence. After iterating over the alignments the script
# will go through this dictionary and create the final ancestor sequence by
# filling up missing positions with gaps.
if len(anc_pos_seq_dict) == 0:
	sys.exit(f"ERROR: No ancestral sequences found in {args.input}!\n"
			f"There has been a configuration/software error or there are "
			f"no blocks for which an ancestral sequence could be obtained.\n"
			f"Please check that the filter includes the ancestral_sequence, "
			f"or whether any were marked in the first place.")
print(f"# nr of sequences found: {len(anc_pos_seq_dict)}, starting at pos "
	f"{next(iter(anc_pos_seq_dict.keys()))}")
# Removes the key for which the length of the ancestral seq is not the same
# as in its annotation.
for key, value in list(anc_pos_seq_dict.items()):
	# if 'Removed' in value:
	if str(value) == 'Removed':
		del anc_pos_seq_dict[key]

# The genomic positions are stored and then sorted, such as the smallest comes first. 
# It will then iterate over this list and fill up a SeqRecord. 
# Before filling of that seqrecord it will insert as many "-" as current_loc - (previous_loc+previous_size)
pos_list = anc_pos_seq_dict.keys()
pos_list_s = sorted(pos_list)

pregaps = '-' * list(pos_list_s)[0]
# pregaps = '-' * (pos_list_s[0]) #ADDED AS MAF IS 0-BASED SO START POS is +1 and pregaps is pos

ancestor_record = SeqIO.SeqRecord(id='',seq=Seq(''))
ancestor_record.seq = ancestor_record.seq + Seq(pregaps)

# Adds the first position
ancestor_record.seq = ancestor_record.seq + anc_pos_seq_dict[list(pos_list_s)[0]].seq

# Enumerate through pos_lis_s starting at second position, thus the index i is always pointing to the previous position
for i,position in enumerate(list(pos_list_s)[1:]):
	# The start of the previous position needs to be bigger than the ancestor_record
	
	# The start of the previous position needs to be bigger than the ancestor_record
	if position+1 < len(ancestor_record):
		print('The start position of ' + str(anc_pos_seq_dict[position].annotations['start']) + ' is located in the previously added sub sequence that is already appended to the Ancestral sequence!')
		prev_position = list(pos_list_s)[i]
		prev_seq = anc_pos_seq_dict[prev_position].seq
		curr_seq = anc_pos_seq_dict[position].seq

		# Find the overlapping part between the previous and current sequences
		if len(ancestor_record) - position > len(curr_seq):
			overlap_len = len(curr_seq)
			print("sequence within prev " +str(overlap_len))
		else:
			overlap_len = len(ancestor_record)-position #prev_position + (len(ancestor_record)-position) - position #
			overlap = str(ancestor_record.seq).endswith(str(curr_seq[:overlap_len]))
			print("#"+str(overlap_len) + " chars, same chars?: " + str(overlap))
			print("Overlap: "+ str(overlap))

		# Append the non-overlapping part of the current sequence to the existing ancestral record
		if overlap:
			non_overlap_seq = curr_seq[overlap_len:]
			ancestor_record.seq += non_overlap_seq
		else:
			overlap_seq = curr_seq[:overlap_len]
			non_overlap_seq = curr_seq[overlap_len:]
			current_fasta = str(ancestor_record.seq)

			curr_start = curr_seq[:overlap_len]
			prev_end = current_fasta[-overlap_len:]
			prev_end = current_fasta[position:position+overlap_len]
			# why are there 2?

			non_matching_indices = []
			for i in range(len(curr_start)):
				if curr_start[i] != prev_end[i]:
					non_matching_indices.append(i)
			print("non matching indices: " + str(non_matching_indices))
			print("Replaced with gaps: " + str(len(non_matching_indices)))
			working_pos = position #len(current_fasta[:prev_position+len(prev_seq)-len(prev_end)])

			for i in non_matching_indices:
				current_fasta = current_fasta[:working_pos+int(i)] + "-" + current_fasta[working_pos+int(i)+1:]
				current_fasta += non_overlap_seq
			ancestor_record.seq = current_fasta
	else:
		# Adds gaps, marking the transition from the previous to the current ancestor sequence
		pregaps = '-'*int(position-(len(ancestor_record))) # SAME AS list(pos_list_s)[i]+anc_pos_seq_dict[list(pos_list_s)[i]].annotations['size'])), SHOULD IT BE -1?
		print("number of pregaps: "+ str(len(pregaps)))
		ancestor_record.seq = ancestor_record.seq + Seq(pregaps)

		# The processed ancestor sequence is appended to the gaps (gaps within the seq are automatically removed)
		ancestor_record.seq = ancestor_record.seq + anc_pos_seq_dict[position].seq

# Adds gaps until chromosome size is reached
ancestor_record.seq = ancestor_record.seq + Seq('-'*(int(chr_lengths[chromosome]) - len(ancestor_record.seq)))

print(f"## Annotating and writing ancestral sequence to file: {args.output}")
# Adds annotations
src_size = anc_pos_seq_dict[position].annotations['srcSize']
ancestor_record.annotations = {'start': 0, 'srcSize': src_size, 'strand': 1, 'size': len(ancestor_record)}
ancestor_record.name = anc_pos_seq_dict[position].name
ancestor_record.id = anc_pos_seq_dict[position].id

ancestor_record.description = f'start: 0, srcSize: {src_size} strand: 1, size: {len(ancestor_record)} '

# Writes ancestor sequence into file
out_f = str(args.output)
with gzip.open(out_f, "wt") if out_f.endswith(".gz") else open(out_f, "w") as output_handle:
	SeqIO.write(ancestor_record, output_handle, "fasta")
