#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys 
from Bio import AlignIO
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


maffile = sys.argv[1]
fastafile = sys.argv[2]
formattedfile = sys.argv[3]
linearizedfile = sys.argv[4]

if(len(sys.argv) != 5):
        sys.exit("usage: format_pairwise_alignments.py MAF_MSA INTERMEDIATE_FASTA INTERMEDIATE_FORMATTED_FASTA OUTPUT_FASTA")

input_handle = open(maffile, "r")
output_handle = open(fastafile, "w")

alignments = AlignIO.parse(input_handle, "maf")
AlignIO.write(alignments, output_handle, "fasta")

output_handle.close()
input_handle.close()

ifile = open(fastafile, 'r')
ofile = open(formattedfile, 'w')

for line in ifile:
	if line.startswith('>'):
		ofile.write('\n' + line)
	else:
		if len(line.strip()) != 60:
			line = line.strip() + '-'*(60-len(line.strip()))
			ofile.write(line.strip())
		else:
			ofile.write(line.strip())

ifile.close()
ofile.close()

species = defaultdict(list)
ifile = open(formattedfile, 'r')
ofile = open(linearizedfile, 'w')

for record in SeqIO.parse(formattedfile, "fasta"):
    record.id = record.id.split('.')[0]
    if record.id not in species:
        species[record.id] = []
    species[record.id].append(record)

ik = dict()
for i, k in enumerate(species):
    ik[i] = k   # dictionary key_of_index

for i in range(0, len(ik)):
    for j in range(0, len(species[ik[0]])):
        if j == 0:
            if i == 0:
                ofile.write('>' + ik[i] + '\n')
            else: 
                ofile.write('\n>' + ik[i] + '\n')
        else:
            ofile.write(str(species[ik[i]][j].seq))


