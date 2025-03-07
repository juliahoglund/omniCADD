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
import copy
from argparse import ArgumentParser
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import logging

# Set up logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

parser = ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input",
	help="maf alignment file from which the ancestral sequence is to be extracted. It is expected that it has been filtered and sorted.", 
	required=True)
parser.add_argument("-o", "--output", 
	help="Ancestral sequence output file, fasta format", 
	required=True)
parser.add_argument("-a", "--ancestor", 
	help="Name/label of the ancestral sequence", 
	required=True)
parser.add_argument("-n", "--interest", 
	help="Name/label of the species of interest",
	required=True)
parser.add_argument("-r", "--reference", 
	help="Name of the genomewide reference species fasta sequence file",
	required=True)

args = parser.parse_args()

def ancestor_finder(alignment_instance, anc_identifer):
    """
    Checks if the ancestor sequence of interest is available in alignment block and returns seq, otherwise empty string.
    :param alignment_instance: Alignment block
    :param anc_identifer: Identifier for the ancestor sequence
    :return: Ancestor sequence or None
    """
    logging.debug(f"Searching for ancestor in alignment block with identifier: {anc_identifer}")
    for entrance in alignment_instance:
        logging.debug(f"Checking sequence: {entrance.id}")
        if entrance.id.startswith(anc_identifer):
            logging.debug(f"Found ancestor sequence: {entrance.id}")
            return entrance
    logging.debug("Ancestor sequence not found in this alignment block.")
    return None

def sequence_processing(g_seq, ancestor_seq):
    """
    Removes gaps that are present in both the sequence of the given species and its ancestor.
    :param g_seq: Sequence of the given species
    :param ancestor_seq: Ancestor sequence
    :return: Processed ancestor sequence or 'Removed' if length mismatch
    """
    logging.debug(f"Processing sequences: {g_seq.id}, {ancestor_seq.id}")
    anc_out_seq = SeqIO.SeqRecord(id=ancestor_seq.id, seq=Seq(''), name=ancestor_seq.name, annotations={'start': ancestor_seq.annotations['start'], 'srcSize': ancestor_seq.annotations['srcSize'], 'strand': ancestor_seq.annotations['strand'], 'size': ancestor_seq.annotations['size']})
    for g_char, anc_char in zip(g_seq, ancestor_seq):
        if g_char != '-':
            anc_out_seq.seq = anc_out_seq.seq + Seq(anc_char)
    if len(anc_out_seq) != anc_out_seq.annotations['size']:
        logging.warning(f"Removed sequence {ancestor_seq.id} due to length mismatch.")
        return 'Removed'
    else:
        return anc_out_seq

def get_alignment_generator(input_f):
    """
    Reads input alignment using AlignIO. Supports either a .maf file or .maf.gz if compressed.
    :param input_f: Path to input file, maf or maf.gz format.
    :return: AlignIO parser on the input file
    """
    logging.debug(f"Reading input file: {input_f}")
    if not os.path.isfile(input_f):
        logging.error(f"File not found: {input_f}")
        sys.exit(f"Program stops prematurely\n Parsed Path does not show to a file. The parsed file is following: {input_f}")
    handle = gzip.open(input_f, "rt") if str(input_f).endswith(".gz") else open(input_f, "r")
    return AlignIO.parse(handle, "maf")

def get_ancestral_sequences(alignment_gen, ancestral_n, interest_n):
    """
    Get dict of ancestral sequences with key start_pos.
    :param alignment_gen: Iterator providing SeqIO records for each alignment block
    :param ancestral_n: Name of ancestral seq
    :param interest_n: Name of the seq of interest
    :return: Dict, key: start_pos, value: ancestral sequence
    """
    logging.debug(f"Extracting ancestral sequences for: {ancestral_n} and {interest_n}")
    anc_dict = {}
    for alignment in alignment_gen:
        sp_seq = alignment[0]
        anc_seq = ancestor_finder(alignment, ancestral_n)
        if (not isinstance(anc_seq, SeqIO.SeqRecord)) or (str(interest_n) + '.' not in sp_seq.id):
            if (isinstance(anc_seq, SeqIO.SeqRecord)) and (str(interest_n) + '.' not in sp_seq.id):
                logging.warning(f"Alignment block contains ancestor sequence but not {interest_n} at first position. Skipping block.")
                continue
            else:
                continue
        elif not isinstance(sequence_processing(sp_seq, anc_seq), SeqIO.SeqRecord):
            continue
        else:
            anc_seq_processed = sequence_processing(sp_seq, anc_seq)
            ref_start = sp_seq.annotations['start']
            if ref_start in anc_dict:
                curr_seq = anc_seq_processed.seq
                other_seq = anc_dict[ref_start].seq
                updated_seq = ''
                seq_length = min(len(curr_seq), len(other_seq))
                for i in range(seq_length):
                    if curr_seq[i] != other_seq[i]:
                        updated_seq += '-'
                    else:
                        updated_seq += curr_seq[i]
                if len(curr_seq) > len(other_seq):
                    updated_seq += curr_seq[seq_length:]
                elif len(other_seq) > len(curr_seq):
                    updated_seq += other_seq[seq_length:]
                updated_seq_record = copy.deepcopy(anc_dict[ref_start])
                updated_seq_record.seq = Seq(updated_seq)
                anc_dict[ref_start] = updated_seq_record
            else:
                anc_dict[ref_start] = anc_seq_processed
    return anc_dict

def get_chr(input_n):
    """
    Get chromosome from input filename.
    :param input_n: Input filename
    :return: Chromosome identifier
    """
    if 'chr' in input_n:
        chr_num = input_n.split('chr')
        chromosome = chr_num[1].split('.')[0]
    else:
        chromosome = 'N'
    return chromosome

def get_chr_lengths(reference_fai):
    """
    Get correct source chromosome size for reference.
    :param reference_fai: Reference fasta index file
    :return: Dict of chromosome lengths
    """
    chr_lengths = defaultdict(dict)
    with open(reference_fai) as f:
        lines = f.read().splitlines()
    lines = [x for x in lines if len(x.split('\t')[0]) < 6]
    for line in lines:
        chrom = line.split()[0]
        length = line.split()[1]
        chr_lengths[chrom] = length
    return chr_lengths


# Read alignment using AlignIO, extract ancestral sequences and put them in
# dict start_pos=sequence. After iterating over the alignments the script
# will go through this dictionary and create the final ancestor sequence by
# filling up missing positions with gaps
print("## Reading MAF file")
alignment_generator = get_alignment_generator(args.input)
anc_pos_seq_dict = get_ancestral_sequences(alignment_generator, args.ancestor, args.interest)
print("## Finished loading ancestral sequences from alignment")

chromosome = get_chr(args.input)
chr_lengths = get_chr_lengths(args.reference)


# Read alignment using AlignIO, extract ancestral sequences and put them in
# dict start_pos=sequence. After iterating over the alignments the script
# will go through this dictionary and create the final ancestor sequence by
# filling up missing positions with gaps.
if len(anc_pos_seq_dict) == 0:
    logging.error(f"No ancestral sequences found in {args.input}!")
    sys.exit(f"ERROR: No ancestral sequences found in {args.input}!\n"
             f"There has been a configuration/software error or there are "
             f"no blocks for which an ancestral sequence could be obtained.\n"
             f"Please check that the filter includes the ancestral_sequence, "
             f"or whether any were marked in the first place.")

print(f"# nr of sequences found: {len(anc_pos_seq_dict)}, starting at pos {next(iter(anc_pos_seq_dict.keys()))}")

# Removes the key for which the length of the ancestral seq is not the same
# as in its annotation.
for key, value in list(anc_pos_seq_dict.items()):
    if str(value) == 'Removed':
        del anc_pos_seq_dict[key]

# The genomic positions are stored and then sorted, such as the smallest comes first. 
# It will then iterate over this list and fill up a SeqRecord. 
# Before filling of that seqrecord it will insert as many "-" as current_loc - (previous_loc+previous_size)
pos_list = anc_pos_seq_dict.keys()
pos_list_s = sorted(pos_list)

pregaps = '-' * list(pos_list_s)[0]
ancestor_record = SeqIO.SeqRecord(id='', seq=Seq(''))
ancestor_record.seq = ancestor_record.seq + Seq(pregaps)

# Adds the first position
ancestor_record.seq = ancestor_record.seq + anc_pos_seq_dict[list(pos_list_s)[0]].seq

for i, position in enumerate(list(pos_list_s)[1:]):
    if position + 1 < len(ancestor_record):
        logging.warning(f"The start position of {anc_pos_seq_dict[position].annotations['start']} is located in the previously added sub sequence that is already appended to the ancestral sequence!")
        prev_position = list(pos_list_s)[i]
        prev_seq = anc_pos_seq_dict[prev_position].seq
        curr_seq = anc_pos_seq_dict[position].seq

        if len(ancestor_record) - position > len(curr_seq):
            overlap_len = len(curr_seq)
        else:
            overlap_len = len(ancestor_record) - position
            overlap = str(ancestor_record.seq).endswith(str(curr_seq[:overlap_len]))

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
            prev_end = current_fasta[position:position + overlap_len]

            non_matching_indices = []
            for i in range(len(curr_start)):
                if curr_start[i] != prev_end[i]:
                    non_matching_indices.append(i)

            working_pos = position

            for i in non_matching_indices:
                current_fasta = current_fasta[:working_pos + int(i)] + "-" + current_fasta[working_pos + int(i) + 1:]
                current_fasta += non_overlap_seq
            ancestor_record.seq = current_fasta
    else:
        pregaps = '-' * int(position - len(ancestor_record))
        ancestor_record.seq = ancestor_record.seq + Seq(pregaps)

        # The processed ancestor sequence is appended to the gaps (gaps within the seq are automatically removed)
        ancestor_record.seq = ancestor_record.seq + anc_pos_seq_dict[position].seq

# Adds gaps until chromosome size is reached
ancestor_record.seq = ancestor_record.seq + Seq('-'*(int(chr_lengths[chromosome]) - len(ancestor_record.seq)))

print(f"## Annotating and writing ancestral sequence to file: {args.output}")
src_size = anc_pos_seq_dict[position].annotations['srcSize']
ancestor_record.annotations = {'start': 0, 'srcSize': src_size, 'strand': 1, 'size': len(ancestor_record)}
ancestor_record.name = anc_pos_seq_dict[position].name
ancestor_record.id = anc_pos_seq_dict[position].id

ancestor_record.description = f'start: 0, srcSize: {src_size} strand: 1, size: {len(ancestor_record)} '

out_f = str(args.output)
with gzip.open(out_f, "wt") if out_f.endswith(".gz") else open(out_f, "w") as output_handle:
    SeqIO.write(ancestor_record, output_handle, "fasta")
