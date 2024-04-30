#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys 
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
   
species = defaultdict(list)
entries = defaultdict(list)

ffile = sys.argv[1] # fasta file to be formatted
output=sys.argv[2] # name of formatted file
indexfile=sys.argv[3] # name of index file with genomic positions

if(len(sys.argv) != 4):
        sys.exit("usage: format_alignments.py CONVERTED_MULTI_FASTA_FILE OUTPUT_FILE OUTPUT_INDEXFILE")

# retrieve all species present in any block to get maximum number of aligned species
for record in SeqIO.parse(ffile, "fasta"):
    record.id = record.id.split('.')[0]
    if record.id not in species:
        species[record.id] = []
    entries[record.id].append(len(record))
    species[record.id].append(record)

# prepare additional dictionaries
# ki = dict()
ik = dict()
for i, k in enumerate(entries):
    # ki[k] = i  # dictionary index_of_key
    ik[i] = k   # dictionary key_of_index

# add gaps fo the length of the block, where that specific species is lacking
# make sure all species with lacking blocks gets a sequence of gaps corresponding to the
# size of the block
offset = 1
for i in range(0, len(ik)-1):
    next_species = ik[i + offset]
    for j in range(0, len(entries[ik[0]])-1):
        if entries[ik[0]][j] != entries[next_species][j]:
            entries[next_species].insert(j, entries[ik[0]][j])
            gaps = '-'*(entries[ik[0]][j])
            species[next_species].insert(j, SeqRecord(Seq(gaps)))
        elif j == len(entries[next_species])-1:
            entries[next_species].append(entries[ik[0]][j+1])
            gaps = '-'*(entries[ik[0]][j+1])
            species[next_species].append(SeqRecord(Seq(gaps)))
        else:
            continue

# as of now hardcoded to trailing unwanted '=' characters'
for i in range(0, len(ik)-1):
    for j in range(0, len(entries[ik[0]])):
        if '=' in species[ik[i]][j].seq:
            species[ik[i]][j].seq = species[ik[i]][j].seq[:-1]

# write linearized fasta sequence
counter = defaultdict(int)
fasta_seq = open(output, 'w')
for i in range(0, len(ik)):
    for j in range(0, len(entries[ik[0]])):
        if j == 0:
            if i == 0:
                fasta_seq.write('>' + ik[i] + '\n')
            else: 
                fasta_seq.write('\n>' + ik[i] + '\n')
        else:
            fasta_seq.write(str(species[ik[i]][j].seq))

# make indexfile of bp positions in reference species
index_out = open(indexfile, 'w')
for j in range(1, len(species[ik[0]])):
    index_out.write(
        'start: ' + str(species[ik[0]][j].description.split(" ")[1]) + 
        ', size: ' + str(species[ik[0]][j].description.split(" ")[2] + '\n')
        )

####
# check if: 
# sequence lengths end up unequal 
# for i in range(0, len(ik)-1):
#     for index, (first, second) in enumerate(zip(entries['sus_scrofa'], entries[ik[i]]), start=1):
#             if first != second:
#                 print('sus scrofa versus', ik[i])
#                 print(index, first, second)
