#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict

def main():
    parser = ArgumentParser(description="Format alignments")
    parser.add_argument("ffile", help="Converted multi-fasta file")
    parser.add_argument("output", help="Name of formatted file")
    parser.add_argument("indexfile", help="Name of index file with genomic positions")
    args = parser.parse_args()

    species = defaultdict(list)
    entries = defaultdict(list)

    # retrieve all species present in any block to get maximum number of aligned species
    for record in SeqIO.parse(args.ffile, "fasta"):
        record.id = record.id.split('.')[0]
        if record.id not in species:
            species[record.id] = []
        entries[record.id].append(len(record))
        species[record.id].append(record)

    # prepare additional dictionaries
    ik = dict()
    for i, k in enumerate(entries):
        ik[i] = k   # dictionary key_of_index

    # add gaps for the length of the block, where that specific species is lacking
    offset = 1
    for i in range(0, len(ik)-1):
        next_species = ik[i + offset]
        for j in range(0, len(entries[ik[0]])-1):
            if entries[ik[0]][j] != entries[next_species][j]:
                entries[next_species].insert(j, entries[ik[0]][j])
                gaps = '-' * entries[ik[0]][j]
                species[next_species].insert(j, SeqRecord(Seq(gaps)))
            elif j == len(entries[next_species])-1:
                entries[next_species].append(entries[ik[0]][j+1])
                gaps = '-' * entries[ik[0]][j+1]
                species[next_species].append(SeqRecord(Seq(gaps)))
            else:
                continue

    # as of now hardcoded to trailing unwanted '=' characters
    for i in range(0, len(ik)-1):
        for j in range(0, len(entries[ik[0]])):
            if '=' in species[ik[i]][j].seq:
                species[ik[i]][j].seq = species[ik[i]][j].seq.rstrip('=')

    # write linearized fasta sequence
    with open(args.output, 'w') as fasta_seq:
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
    with open(args.indexfile, 'w') as index_out:
        for j in range(1, len(species[ik[0]])):
            index_out.write(
                'start: ' + str(species[ik[0]][j].description.split(" ")[1]) +
                ', size: ' + str(species[ik[0]][j].description.split(" ")[2]) + '\n'
            )

if __name__ == "__main__":
    main()
