#!/usr/bin/python3

# Author: Tom van der Valk
# Edit: julia höglund, 2025-06-12
# Modification: This script has been modified to add a column with base pair positions

from sys import argv
from itertools import islice
import operator
import gzip

def parse_gerp(fasta_file, gerp_file, species_name):

    seq = ""
    neutral_score = []
    RS_score = []
    outputfile = open(gerp_file + ".parsed", "w")

    adder = False
    with open(fasta_file) as f1:
        for line in f1:
            if line.startswith(">"):
                if line.startswith(">" + species_name):
                    adder = True
                else:
                    adder = False
            if adder and not line.startswith(">"):
                seq += line.strip()

    N_count = seq.count("N")

    with open(gerp_file) as f2:
        for line in f2:
            splitted = line.strip().split("\t")
            neutral_score += [splitted[0]]
            RS_score += [splitted[1]]

    position = 0
    gerp_position = 0
    for i in seq:
        position += 1
        if i.upper() == "N":
            outputfile.write(f"{position}\t0\t0\n")
            gerp_position += 1
        else:
            idx = position - gerp_position - 1
            if idx < len(neutral_score):
                outputfile.write(f"{position}\t{neutral_score[idx]}\t{RS_score[idx]}\n")

    outputfile.close()
    print('file ' + str(fasta_file) + ' done.')

if __name__ == "__main__":
    fasta_file = argv[1]
    gerp_file = argv[2]
    species_name = argv[3]
    parse_gerp(fasta_file, gerp_file, species_name)