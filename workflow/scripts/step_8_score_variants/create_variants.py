#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Christian Geoss
:Contact: c.gross@tudelft.nl
:Date: 06.05.2019

:Edited by: Job van Schipstal
:Date: 24-10-2023
:Usage: see create_variants.py --help

This script iterates over the reference sequence and generates all possible
SNPs for the chromosome of interest. These are stored in .vcf files of a
user specified size at the desired location.
"""

import re
import sys
from argparse import ArgumentParser

from Bio import SeqIO

parser = ArgumentParser(__doc__)
parser.add_argument("-o", "--out-path", type=str,
                    help="Path to output folder, "
                         "(default: the working directory)", default="")
parser.add_argument("-s", "--size", type=int,
                    help="Number of sites for which variants are created. "
                         "The total number of lines in the output is 3 times "
                         "the specified number, for the 3 alternative alleles "
                         "per site. (default: 1000000)", default=1000000)
parser.add_argument("-r", "--reference", type=str,
                    help="Path to reference, fasta format", required=True)
parser.add_argument("-p", "--start-pos", type=int,
                    help="The position at which to start generating variants, "
                         "0-based (default: 0)", default=0)
parser.add_argument("-c", "--chrom", type=str,
                    help="The chromosome for which to generate variants, "
                         "the name of the chromosome must be in the reference "
                         "fasta file (default: 1)", default="1")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase verbosity, prints to stdout")
args = parser.parse_args()

NUCLEOTIDES = ["A", "C", "G", "T"]


def get_chrom_from_file(chromosome, fasta_file, vb=False):
    """
    Attempts to find the desired sequence in a fasta file,
    exits if unsuccessful.

    Chromosome must be part of the name, either next to white space,
    at the end/start of the line or prefixed by chr or chromosome.

    :param chromosome: str, name of chromosome of interest
    :param fasta_file: file to read from
    :param vb: optional boolean, should print verbose output, default False
    :return: SeqIO record for chromosome of interest
    """
    pattern = r"(^|chr|chromosome)\W*" + chromosome + r"($|\W)"
    ref_fasta = SeqIO.parse(fasta_file, "fasta")
    for rec in ref_fasta:
        return rec  # TODO CHANGE, EITHER CHECK BETTER OR NOT
        if vb:
            print(f"Found: ID = {rec.id}, length {len(rec.seq)}, "
                  f"with {len(rec.features)} features")
        if re.match(pattern, rec.name):
            print(f"Chromosome found:\nName: {rec.name}\n"
                  f"Length {len(rec.seq)}")
            return rec
    sys.exit(f"Chr {chromosome} was not found in reference fasta file:\n"
             f"{fasta_file}")


def open_vcf_out(path, chrom, part, vb=False):
    """
    Opens a part file in the desired folder and writes VCF header.

    :param path: str, prefix for the file, path to save output in.
    :param chrom: str, chromosome, included in name
    :param part: int, current part, included in name
    :param vb: optional boolean, print debug information
    :return: openfile instance, vcf file with header
    """
    if not path.endswith("/"):
        path += "/"
    file = f"{path}/{part}.vcf"
    if vb:
        print(f"Opening file: {file}")
    file_h = open(file, "w")
    file_h.write("##fileformat=VCFv4.1\n")
    file_h.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    return file_h


def jump_to_pos(enumerate_iter, stop_index, vb=False):
    """
    Iterates the enumerator to the desired position, does nothing for index <=0

    :param enumerate_iter: instance of enumarate to jump in.
    :param stop_index: int, index to jump to
    :param vb: optional boolean, should print debug info to stdout
    :return: input enumarate instance, jummped to desired position
    """
    if stop_index <= 0:
        return enumerate_iter
    if vb:
        print("Jumping to start site on identified chromosome")
    for i, j in enumerate_iter:
        if i == stop_index - 1:
            return enumerate_iter
        elif (i % 10000000 == 1) and vb:
            print(f"Current Position {i}")


record = get_chrom_from_file(args.chrom, args.reference, args.verbose)

if args.start_pos >= len(record.seq):
    sys.exit(f"The specified start location is outside of the length of the "
             f"specified start chromosome.\n Specified location is: "
             f"{args.coordinates}\n Chromosome length is "
             f"{record.id}:{len(record.seq)}")

# Open the first outfile, new ones will be opened when full
outfile_iter = 1
outfile = open_vcf_out(args.out_path, args.chrom, outfile_iter, args.verbose)

# Get sequence enumerator, starting at the desired position
iter_enumerate = jump_to_pos(enumerate(record.seq),
                             args.start_pos, args.verbose)

# CAUTION: site_index is 0-based but for the vcf file it has to be 1-based
remaining_size = args.size
for site_index, ref_allele in iter_enumerate:
    # Ignore gaps or uncertain nt's
    ref_allele = ref_allele.upper()
    if ref_allele not in NUCLEOTIDES:
        continue

    for alt_allele in NUCLEOTIDES:
        if ref_allele != alt_allele:
            outfile.write(f"{args.chrom}\t{site_index + 1}\t"
                          f".\t{ref_allele}\t{alt_allele}\t.\t.\t.\n")
    remaining_size -= 1

    # Open new file if chunk size is reached
    if remaining_size <= 0:
        outfile.close()
        outfile_iter += 1
        remaining_size = args.size
        outfile = open_vcf_out(args.out_path, args.chrom,
                               outfile_iter, args.verbose)

outfile.close()
print(f"# Finished generating variants.\n"
      f"{(outfile_iter * args.size - remaining_size) * 3} variants.\n"
      f"{outfile_iter} parts.")
