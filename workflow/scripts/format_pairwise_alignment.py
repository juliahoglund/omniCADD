#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from Bio import AlignIO, SeqIO
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def convert_maf_to_fasta(maffile: str, fastafile: str) -> None:
    """
    Convert MAF file to FASTA format.
    """
    with open(maffile, "r") as input_handle, open(fastafile, "w") as output_handle:
        alignments = AlignIO.parse(input_handle, "maf")
        AlignIO.write(alignments, output_handle, "fasta")

def format_fasta(fastafile: str, formattedfile: str) -> None:
    """
    Format FASTA file to ensure each sequence line is 60 characters long.
    """
    with open(fastafile, 'r') as ifile, open(formattedfile, 'w') as ofile:
        for line in ifile:
            if line.startswith('>'):
                ofile.write('\n' + line.strip() + '\n')
            else:
                line = line.strip()
                ofile.write(line + '-' * (60 - len(line)) + '\n')

def linearize_fasta(formattedfile: str, linearizedfile: str) -> None:
    """
    Linearize FASTA file by concatenating sequences of the same species.
    """
    species = defaultdict(list)
    with open(formattedfile, 'r') as ifile:
        for record in SeqIO.parse(ifile, "fasta"):
            record.id = record.id.split('.')[0]
            species[record.id].append(record)

    with open(linearizedfile, 'w') as ofile:
        ik = {i: k for i, k in enumerate(species)}
        for i in range(len(ik)):
            for j in range(len(species[ik[0]])):
                if j == 0:
                    if i == 0:
                        ofile.write('>' + ik[i] + '\n')
                    else:
                        ofile.write('\n>' + ik[i] + '\n')
                else:
                    ofile.write(str(species[ik[i]][j].seq))

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit("usage: format_pairwise_alignments.py MAF_MSA INTERMEDIATE_FASTA INTERMEDIATE_FORMATTED_FASTA OUTPUT_FASTA")

    maffile, fastafile, formattedfile, linearizedfile = sys.argv[1:5]

    try:
        convert_maf_to_fasta(maffile, fastafile)
        format_fasta(fastafile, formattedfile)
        linearize_fasta(formattedfile, linearizedfile)
    except Exception as e:
        print(f"Error processing files: {e}")
        sys.exit(1)


