#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 06-08.2018

This files goes through the ancestral sequence 
and compares it with the reference genome of interest
and the known population variants to identfy derived variants.

The identification is based on five criteria:
Cases		:1		2		3		4		5
Ancestor	:C		C		C		C		C
Reference	:C		A>0.9	A		A>0.9	A
Alternative	:A>0.9	G		G>0.9	C		A

Case 1: The reference genome and ancestor are the same at position n, 
		but in the population data, there is a nearly fixed derived variant with >0.9
Case 2: The reference genome and ancestor are different at position n, 
		and the reference allele has a frequency >0.9
Case 3: The reference genome and ancestor are different at position n, 
		but there is a third allele in the population with a frequency >0.9
Case 4: The reference genome and ancestor are different, but the alternative allele
		in the population data is not the one in the reference genome, i.e. is most 
		likely ancestral. Here the population reference allele at >0.9 is what is desired
Case 5: The ancestral and reference allele are different at position n

Case X: The reference genome and ancestor are different, but there is no data on that position
		in the population file. If a true variant and true missingness, it can be assumed to be 100% fixed
		in the reference species / population, and hence be monomorphic in the reference.
		not yet implemented

:Edited by: Seyan Hu
:Date: 19-10-2022
:Extension and modification: Julia Höglund
:Date: 2023-08-08
:Edited by: Job van Schipstal
:Date: 27-10-2023
:Usage: See derive_variants.py --help
Changed to better follow PEP8 style and switched to argparse for input.
"""

import gzip
import sys
from argparse import ArgumentParser
import pysam

parser = ArgumentParser(description=__doc__)
parser.add_argument('-o', '--output',
	help='Vcf output file prefix and path, both low and high quality variants will be generated '
	'(default: derived_variants_)',
	default="derived_variants_")
parser.add_argument('-c', '--chrom',
	help='chromosome of interest',
	required=True)
parser.add_argument('-a', '--ancestor',
	help='The ancestral sequence file for the chromosome',
	required=True)
parser.add_argument('-r', '--reference',
	help='The reference genome file for the chromosome',
	required=True)
parser.add_argument('-v', '--variants',
	help='The population variant file for the chormosome',
	required=True)
parser.add_argument('-s', '--start', 
	help='start position of the region to generate variants for',
	default=0, 
	type=int)

# Check if the nt are written in upper (unmasked) or lower case (masked), 
# depending on this, the script writes it in one of the two output files.
def write_line(chrom, pos, reference, an, output):
	if (reference in 'ACGTacgt') and (an in 'ACGTacgt'):
		if an.islower():
			output.write(f"{chrom}\t{pos + 1}\t.\t{reference}\t{an.upper()}\t.\t.\t.\n")
		else:
			output.write(f"{chrom}\t{pos + 1}\t.\t{reference}\t{an}\t.\t.\t.\n")

def main(args):
	# Define input files for the Ancestral sequence and the reference.
	ancestor_fasta = pysam.Fastafile(args.ancestor)
	ref_fasta = pysam.Fastafile(args.reference)

	# Open population frequency files
	# Note: frequency files are not vcf files
	input_snps = gzip.open(args.variants, "rt") \
		if args.variants.endswith('.gz') else open(args.variants, "r")

	# Check if the ancestor seq is of same length as the reference
	if ancestor_fasta.nreferences != 1:
		sys.exit('There are more than one (1) record in the fasta file.')

	ancestor_record = ancestor_fasta.references[0]

	if f'.{args.chrom}' in ancestor_record or f'chr{args.chrom}' in ancestor_record:
		anc_length = ancestor_fasta.get_reference_length(ancestor_fasta.references[0])
	else:
		sys.exit(
			'The requested chromosome cannot be found in the record of '
			'the ancestor fasta\n%s\n%s' % (f'.{args.chrom}:', ancestor_record))

	if anc_length != ref_fasta.get_reference_length(ref_fasta.references[0]):
		sys.exit(
			'Ancestor fasta and ref fasta does not have the same length, '
			'Chromosome: %s' % (f'{args.chrom}'))

	print(f'Ancestor sequence control: Chr{args.chrom}, sequence is good')

	# Create output files for upper and lower case variants.
	output = open(f"{args.output}.vcf", 'w')
	output.write("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

	# Split the ancestor and reference sequences in chunks of 350000.
	# Done so that the script reads this in chunks.
	index_outer = (args.start, args.start + 350000)
	input_snps.readline()  # Skip header
	snp_line = input_snps.readline()
	snp_pos = int(snp_line.strip().split('\t')[1])

	while index_outer[0] <= anc_length:
		# Fetch the reference and ancestral sequence.
		ref_seq = ref_fasta.fetch(ref_fasta.references[0], index_outer[0], index_outer[1])
		anc_seq = ancestor_fasta.fetch(ancestor_record, index_outer[0], index_outer[1])

		# Iterate over the ancestral sequence position. 
		for i, anc in enumerate(anc_seq):
			
			# Enumerate sets i always to 0, therefore i needs to be summed with the starting position of index_outer
			# Checking if current position in the sequences is equal to the snp position in the frequency files
			if i + index_outer[0] == (snp_pos-1):
				snp_allele = snp_line.strip().split('\t')[5].split(':')[0].upper()
				
				# Checking if Case 1 derived snp is given
				if (anc.upper() == ref_seq[i]) and (anc.upper() != snp_allele) and (anc.upper() != '-') and (anc.upper() != '.'):
					# Case 1 write out
					write_line(args.chrom, i + index_outer[0], ref_seq[i], snp_allele, output)
					try: 
						snp_line = input_snps.readline()
						snp_pos = int(snp_line.strip().split('\t')[1])
					except:
						pass
				
				# Checking if case 3 (or 2) derived snp is given
				elif (anc.upper() != ref_seq[i]) and (anc.upper() != snp_allele) and (anc.upper() != '-'):
					# Case 3 write out
					write_line(args.chrom, i + index_outer[0], ref_seq[i], snp_allele, output)
					try: 
						snp_line = input_snps.readline()
						snp_pos = int(snp_line.strip().split('\t')[1])
					except: 
						pass
				
				else: 
					try:
						snp_line = input_snps.readline()
						snp_pos = int(snp_line.strip().split('\t')[1])
					except: 
						pass
						
			elif anc.upper() != '-' :
				# Checking if case 5 (or 3) derived variant
				if (anc.upper() != ref_seq[i]):
					# Case 5 write out
					write_line(args.chrom, i + index_outer[0], ref_seq[i], anc, output)
					
		index_outer = (index_outer[1], index_outer[1] + 350000)		
				
	output.close()

if __name__ == '__main__':
	main(parser.parse_args())

