#!/usr/bin/env python
# -*- coding: ASCII -*-
"""
:Author: Christian Gross
:Contact: cgross@tudelft.nl
:Date: 01-08-2018

This script takes the vcf file with derived variants and filters out those
who lies on consecutive positions, as those might be interpreted as indels,
i.e. longer than SNPs

:Edited by: Seyan Hu
:Date: 20-10-2022
:Edited by: Job van Schipstal
:Date: 22-9-2023
:Usage: see filter_for_singletons.py --help
Reformatted code in accordance with PEP guidelines and switched to argparse.
"""

# Import dependencies
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument('-i', '--input',
                    help='Input vcf file to filter for SNPs',
                    required=True)
parser.add_argument('--snps',
                    help='Output file path for separated SNPs',
                    required=True)
parser.add_argument('--series',
                    help='Output file for adjacent SNPs',
                    required=True)


def separate_snp(infile, snp_outfile, series_outfile) -> None:
    """
    Opens input and output vcf files, reads variants and determines if there
    are adjacent variants, if they are separated by at least 1 position they
    are written to the snp file, otherwise they are written to the series file.
    :param infile: input vcf file to filter
    :param snp_outfile: filename for singleton SNPs
    :param series_outfile: filename for series SNPs
    :return: None, writes to the output files directly
    """
    # Read vcf file and create two seperate files for snps and indels.
    infile = open(infile, 'r')
    snp_f = open(snp_outfile, 'w')
    series_f = open(series_outfile, 'w')
    for out in [snp_f, series_f]:
        out.write('##fileformat=VCFv4.1\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    # Loop through the lines in vcf file
    previous_pos = 0
    previous_lines = ''
    jump = False
    for lines in infile:

        # Skips header
        if lines.startswith('#'):
            continue

        line = lines.split('\t')
        # If 'previous_line' is an empty string, it appends the current
        # line to the previous line and sets 'previous_pos' to the
        # variants position of the current line. Happens when the loop
        # encounters the first line of information.
        if previous_lines == '':
            previous_lines = lines
            previous_pos = int(line[1])
            continue
        # If the position of the current variant is next to the previous
        # position, it writes the line to the indel file and sets the
        # 'previous_line' to the current line and 'previous_pos' to the
        # variants position of the current line. It also sets 'jump' to
        # True.
        if (int(line[1]) - 1 == previous_pos) or (
                int(line[1]) == previous_pos):
            jump = True
            series_f.write(previous_lines)
            previous_lines = lines
            previous_pos = int(line[1])
            continue
        # If the position of the current variant is not next to the
        # previous position and not the first line, it check whenever
        # 'jump' is still True. If 'jump' is still True, it will be set
        # to False and the line is written to the indel file. The
        # 'previous_line' is set to the current line and 'previous_pos'
        # to the variants position of the current line.

        # Else it writes the line to the snp file and
        # the 'previous_line' is set to the current line
        # and 'previous_pos' to the variants position of the current line.
        if jump:
            jump = False
            series_f.write(previous_lines)
            previous_lines = lines
            previous_pos = int(line[1])
            continue
        snp_f.write(previous_lines)
        previous_lines = lines
        previous_pos = int(line[1])

    if jump:
        series_f.write(previous_lines)
    else:
        snp_f.write(previous_lines)

if __name__ == '__main__':
    args = parser.parse_args()
    separate_snp(args.input, args.snps, args.series)
