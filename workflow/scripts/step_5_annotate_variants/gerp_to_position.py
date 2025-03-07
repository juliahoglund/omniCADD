#!/usr/bin/python3

# Author: Tom van der Valk

from sys import argv, exit
from itertools import islice
import operator
import gzip
import logging
from typing import List

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_gerp(fasta_file: str, gerp_file: str, species_name: str) -> None:
    try:
        sequence = ""
        neutral_scores: List[str] = []
        rs_scores: List[str] = []

        with open(fasta_file) as f1:
            adder = False
            for line in f1:
                if line.startswith(">"):
                    adder = line.startswith(">" + species_name)
                elif adder:
                    sequence += line.strip()

        with open(gerp_file) as f2:
            for line in f2:
                splitted = line.strip().split("\t")
                neutral_scores.append(splitted[0])
                rs_scores.append(splitted[1])

        with open(gerp_file + ".parsed", "w") as outputfile:
            position = 0
            gerp_position = 0
            for nucleotide in sequence:
                position += 1
                if nucleotide.upper() == "N":
                    outputfile.write("0\t0\n")
                    gerp_position += 1
                else:
                    if (position - gerp_position - 1) < len(neutral_scores):
                        outputfile.write(neutral_scores[position - gerp_position - 1] + "\t" + 
                                        rs_scores[position - gerp_position - 1] + "\n")

        logging.info('File %s processed successfully.', fasta_file)

    except Exception as e:
        logging.error('Error processing file: %s', e)
        exit(1)

if __name__ == "__main__":
    if len(argv) != 4:
        logging.error('Usage: python gerp_to_position.py <fasta_file> <gerp_file> <species_name>')
        exit(1)
    fasta_file = argv[1]
    gerp_file = argv[2]
    species_name = argv[3]
    parse_gerp(fasta_file, gerp_file, species_name)