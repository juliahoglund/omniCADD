import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
   
species = defaultdict(list)
entries = defaultdict(list)

ffile = sys.argv[1]
output=sys.argv[2] # folder to save formatted chunks in

if(len(sys.argv) != 3):
        sys.exit("usage: format_alignments.py CONVERTED_MULTI_FASTA_FILE OUTPUT_FILE")

# retrieve all species present in any block to get maximum number of aligned species
for record in SeqIO.parse(ffile, "fasta"):
    if record.id not in species:
        species[record.id].append(SeqRecord(Seq(""), id=record.id, name=record.id))
    entries[record.id].append(len(record))
    species[record.id].append(record)

# prepare additional dictionaries
# ki = dict()
ik = dict()
for i, k in enumerate(entries):
    # ki[k] = i   # dictionary index_of_key
    ik[i] = k   # dictionary key_of_index

# add gaps fo the length of the block, where that specific species is lacking
offset = 1
for i in range(0, len(ik)-1):
    next_species = ik[i + offset]
    for j in range(0, len(entries[ik[0]])-1):
        if j == len(entries[next_species])-1: 
            entries[next_species].append(entries[ik[0]][j+1])
            gaps = '-'*(entries[ik[0]][j+1])
            species[next_species].append(SeqRecord(Seq(gaps)))
        elif entries[ik[0]][j] == entries[next_species][j]:
            continue
        else:
            entries[next_species].insert(j, entries[ik[0]][j])
            gaps = '-'*(entries[ik[0]][j])
            species[next_species].insert(j, SeqRecord(Seq(gaps)))

fasta_seq = open(output, 'w')
for i in range(0, len(ik)-1):
    for j in range(0, len(entries[ik[0]])):
        if j == 0:
            if i == 0:
                fasta_seq.write('>' + ik[i] + '\n')
            else: 
                fasta_seq.write('\n>' + ik[i] + '\n')
        else:
            fasta_seq.write(str(species[ik[i]][j-1].seq))
