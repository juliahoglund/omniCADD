#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:Author: Christian Gross
:Contact: cgross@tudelft.nl
:Date: 01-08-2018

This script takes the vcf file with simulated SNPs and iterates over it,
at each row identifying if the variant is at a position with a known high
quality ancestor or not. Depending on that it is splitting the file.

:Edited by: Seyan Hu
:Date: 4-11-2022
:Edited by: Job van Schipstal
:Date: 2-10-2023
:Usage: see python filter_ancestor_site.py --help
- Modified to accept input, ancestral and output filenames from argparse
  instead of searching in a given folder.
- Do not check chromosome for each variant,
  it is expected for the inputs to match.
- Reformatted to better align with PEP8 style guidelines.
"""

# Import dependencies.
import sys
import gzip
from argparse import ArgumentParser
import pysam

parser = ArgumentParser(description=__doc__)
parser.add_argument('-i', '--input',
	help = 'Input vcf file of variants to filter',
	type = str, 
	required = True)
parser.add_argument('-a', '--ancestor',
	help = 'The ancestral sequence for the chromosome',
	type = str, 
	required = True)
parser.add_argument('-o', '--output',
	help = 'output filename, vcf format.',
	type = str, 
	required = True)


def filter_vcf_by_ancestral(input_vcf, output_vcf, ancestral) -> None:
	"""
	Read input vcf and filter variants by presence of
	an high quality ancestral allele.
	:param input_vcf: str, filename for input file
	:param output_vcf: str, filename for output file
	:param ancestral: str, filename of ancestral sequence fasta file
	:return: None, written to file
	"""
	# Open input/output vcf and ancestral fasta file
	vcf_file = gzip.open(input_vcf, "rt") \
		if input_vcf.endswith('.gz') else open(input_vcf, "r")

	outfile = gzip.open(output_vcf, "wt") \
		if output_vcf.endswith('.gz') else open(output_vcf, "w")
	anc_fasta = pysam.Fastafile(ancestral)
	if anc_fasta.nreferences != 1:
		sys.exit(f"Expected 1 sequence in ancestral fasta, not {anc_fasta.nreferences}")

	reference = anc_fasta.references[0]

	# Iterate over lines in vcf.
	for lines in vcf_file:
		if lines[0] == '#':
			outfile.write(lines)
			continue
		pos = lines.split('\t')[1]
		ancestor = anc_fasta.fetch(reference, int(pos) - 1, int(pos))
		if ancestor in 'ACGT':
			outfile.write(lines)
	outfile.close()
	vcf_file.close()


if __name__ == '__main__':
	args = parser.parse_args()
	filter_vcf_by_ancestral(args.input, args.output, args.ancestor)

