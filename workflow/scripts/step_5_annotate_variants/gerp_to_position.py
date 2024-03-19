#!/usr/bin/python3

# Author: Tom van der Valk

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
    for i in (seq):
        position += 1
        if i.upper() == "N":
            outputfile.write("0\t0" + "\n")
            gerp_position += 1
        else:
            if (position-gerp_position-1) < len(neutral_score):
                outputfile.write(neutral_score[position - gerp_position -1] + "\t" + 
                    RS_score[position - gerp_position -1] + "\n")

    outputfile.close()
    print('file ' + str(filename) + ' done.')

if __name__ == "__main__":
    filename = argv[1]
    species_file = argv[2]
    species_name = argv[3]
    parse_gerp(filename, species_file, species_name)