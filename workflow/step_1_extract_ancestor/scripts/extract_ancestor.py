"""
:Author: Christian Gross
:Contact: cgross@tudelft.nl
:Date: 01.02.2018

This script takes as input the preporcessed maf files,
the ancestor sequence of interest and returns per .maf file (per chromosome) the ancestor sequence in fasta format.

:Edited by: Seyan Hu
:Date: 29-9-2022
:Modified by: Julia Höglund
:Date 22-3-2023

:Usage: python extract_ancestor.py <file> <Ancestor identifier> <scientific name> <fai file>

:Example:
python scripts/extract_ancestor.py output/dir_mS/mS_chr19 chr19 Ancestor_Mouse_Rat mus_musculus Mus_Musculus.fai

Added:
1) It removes ancestral sequences that are not of the same length as stated in their annotation.
2) It checks whenever the start position of the to be added sub sequence is located in the previous added sub sequence.
	If true, this sequence will be skipped.

Needs to be looped, because it only takes one chr.maf file
Solved with this script: 'wrapper_extract_ancestor.py'
"""

# Import dependencies
import sys, os
from optparse import OptionParser
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from sh import gunzip
from collections import defaultdict

# Assign input to variables
path_to_input = sys.argv[1]
ancestor_ident = sys.argv[2]
sc_name = sys.argv[3]
index_file = sys.argv[4]

# Determine chr.
if 'chr' in path_to_input:
	chr_num = path_to_input.split('chr')
	chromosome = chr_num[1].split('.')[0]
else:
	chr_num = 'N'

# A function that checks if the ancestor sequence of interest
# is available in alignment block and returns seq, otherwise empty string
def ancestor_finder(alignment_instance,anc_identifer):
	for i,entrance in enumerate(alignment_instance):
		if entrance.id.startswith(anc_identifer):
			return entrance

	return None

# Function that removes gaps that are present in both the sequence of the given species and its ancestor.
# This is done to find the real genomic location of the ancestor to the reference.
# Otherwise it would result in wrong positions when deriving variants.
def sequence_processing(g_seq,ancestor_seq):

	# Initializing two empty SeqRecord objects which will be filled with the processed given species and ancestor sequences respectively.
	# This way the true length of the ancestor can be identified and it will fill up all the remaining positions with gaps.
	anc_out_seq = SeqIO.SeqRecord(id=ancestor_seq.id,seq=Seq(''),name=ancestor_seq.name,annotations={'start': ancestor_seq.annotations['start'], 'srcSize': ancestor_seq.annotations['srcSize'], 'strand': ancestor_seq.annotations['strand'], 'size': ancestor_seq.annotations['size']})

	for (g_char,anc_char) in zip(g_seq,ancestor_seq):
		if g_char!='-':
			anc_out_seq.seq = anc_out_seq.seq + Seq(anc_char)

	# Checks if the length is the same as the size in annotations
	if len(anc_out_seq)!=anc_out_seq.annotations['size']:
		return 'Removed'
	else:
		return anc_out_seq

# Makes sure that the path shows a file which exists
if (not os.path.isfile(path_to_input)):
	sys.exit('Program stopped prematurely \n Parsed path does not point to a file. The parsed file is following : {0}'.format(path_to_input))

####################

print('Working on file: ' + path_to_input)

# get total chromosome length (in bases)
chr_lengths = defaultdict(dict)

with open(index_file) as f:
    lines = f.read().splitlines()
lines = [ x for x in lines if len(x.split('\t')[0]) < 6 ]
for line in lines:
	chrom = line.strip.split()[0]
	length = line.strip.split()[1]
	chr_lengths[chrom] = length

# Read file
handle = open(path_to_input, "r")
alignment_generator = AlignIO.parse(handle, "maf")

# Init dictionary with sequence start location as key and the sequence as values.
# After iterating over the alignments it will go through this dictionary
# and create the final ancestor sequence by filling up missing positions with gaps.
anc_pos_seq_dict = {}

# Loop through the file and create a dict for the ancestral seq.
for alignment in alignment_generator:
	sp_seq = alignment[0]

	# Search for the ancestor sequence of interest in alignment block and checks if the seq is of same size as in the annotation
	anc_seq = ancestor_finder(alignment,ancestor_ident)

	if (not isinstance(anc_seq,SeqIO.SeqRecord)) or (str(sc_name) +'.' not in sp_seq.id):
		if (isinstance(anc_seq,SeqIO.SeqRecord)) and (str(sc_name) +'.' not in sp_seq.id):
			print('In the preprocessing is something wrong. The alignment block contains an ancestor sequence but not '+ str(path_to_input) +' at first position. Path and sequence details %s \t %s' %(path_to_input,anc_seq.annotations))
			print('This alignment block will be skipped.\n')
			continue
		else:
			continue
	else:
		anc_pos_seq_dict[anc_seq.annotations['start']] = sequence_processing(sp_seq,anc_seq)


# Removes the key for which the length of the ancestral seq is not the same as in its annotation.
for key, value in list(anc_pos_seq_dict.items()):
	#if 'Removed' in value:
	if str(value) == 'Removed':
		del anc_pos_seq_dict[key]


# The genomic positions are stored and then sorted, such as the smallest comes first.
# It will then iterate over this list and fill up a SeqRecord.
# Before filling of that seqrecord it will insert as many "-" as current_loc - (previous_loc+previous_size)
pos_list = anc_pos_seq_dict.keys()
pos_list_s = sorted(pos_list)

pregaps = '-' * list(pos_list_s)[0]

ancestor_record = SeqIO.SeqRecord(id='',seq=Seq(''))
ancestor_record.seq = ancestor_record.seq + Seq(pregaps)

# Adds the first position
ancestor_record.seq = ancestor_record.seq + anc_pos_seq_dict[list(pos_list_s)[0]].seq

# Enumerate at second position, thus the index i is always pointing to the previous position
for i,position in enumerate(list(pos_list_s)[1:]):

	# The start of the previous position needs to be bigger than the last position of the previous sequence
	if anc_pos_seq_dict[position].annotations['start'] < len(ancestor_record):
		print('The start position of ' + str(anc_pos_seq_dict[position].annotations['start']) + ' is located in the previously added sub sequence that is already appended to the Ancestral sequence!')
		continue
	else:
		# Adds gaps, marking the transition from the previous to the current ancestor sequence.
		pregaps = '-'*(anc_pos_seq_dict[position].annotations['start']-(anc_pos_seq_dict[list(pos_list_s)[i]].annotations['start']+anc_pos_seq_dict[list(pos_list_s)[i]].annotations['size']))
		ancestor_record.seq = ancestor_record.seq + Seq(pregaps)

		# The processed ancestor sequence is appended to the gaps (gaps within the seq are automatically removed)
		ancestor_record.seq = ancestor_record.seq + anc_pos_seq_dict[position].seq

# Adds gaps until chromosome size is reached
ancestor_record.seq = ancestor_record.seq + Seq('-'*(int(chr_lengths[chromosome]) - len(ancestor_record.seq)))

# Adds annotations
ancestor_record.annotations = {'chromosome': chromosome, 'start': 0, 'srcSize': int(chr_lengths[chromosome]), 'strand': 1, 'size': len(ancestor_record)}
ancestor_record.name = anc_pos_seq_dict[position].name.split('.')[0]
ancestor_record.id = anc_pos_seq_dict[position].id

ancestor_record.description = 'chromosome: {0}, start: 0, srcSize: {1} strand: 1, size: {2}'.format(chromosome,int(chr_lengths[chromosome]),len(ancestor_record))

# Writes ancestor sequence into file
with open('./' + ancestor_record.id + '_chr' + chromosome + '.fa', "w") as output_handle:
	SeqIO.write(ancestor_record, output_handle, "fasta")
