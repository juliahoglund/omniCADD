#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 06-08.2018

This files goes through the ancestral sequence 
and compares it with the reference genome of interest
and the known population variants to identfy derived variants.

The identification is based on X criteria:
Cases		:1		2		3		4		5
Ancestor	:C		C		C		C		C
Reference	:C		A>0.9	A		A>0.9	A
Alternative	:A>0.9	G		G>0.9	C		A

Case:		:1 		  2 		3 			4
Ancestor 	:C 		  C 		C 			C
Reference 	:G 		  C 		G 			G
Population  :ref/alt  C/G>0.9   G/C>0.1 	no data available at position.
											(likely monomorphic in pop file)

Case 1: The reference genome and ancestor are different at position n
Case 2: The reference genome and ancestor are the same at position n, 
		but in the population data, there is a nearly fixed derived variant with >0.9
Case 3: The reference genome and ancestor are different, but the alternative allele
		in the population data is not the one in the reference genome, i.e. is most 
		likely ancestral. Here the population reference allele at >0.9 is what is desired
Case 4: The reference genome and ancestor are different, but there is no data on that position
		in the population file. If a true variant and true missingness, it can be assumed to be 100% fixed
		in the reference species / population, and hence be monomorphic in the reference.

:Edited by: Seyan Hu
:Date: 19-10-2022
:Usage: python derived_var_gen.py <chr num> <path to ancestor seq> <path to genome> <path to frequency files>

:Example:
python derived_var_gen.py 19 ./ ../../../generate_ancestral_seq/data/genome/ ../../output/dir_freq_f/

Note:
Since this script only performs the generation of derived variants for one chromosome, 
a wrapper is needed to iterate over all chromosomes. 
"""


# Import dependencies
import os,sys
import io
from optparse import OptionParser
import time
import pysam
from pysam import VariantFile
import natsort

# OptionParser for input.
parser = OptionParser()
parser.add_option("-c", "--chromosome", dest = "chromosome", help = "list of chromosomes from which the variants are going to be derived", default = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X")
parser.add_option("-a", "--ancestor", dest = "ancestor", help = "path to the folder with the extracted ancestor", default = "data/extracted_ancestor/")
parser.add_option("-g", "--genome", dest = "genome", help = "path to the reference genome folder", default = "data/genome/")
parser.add_option("-f", "--frequency", dest = "frequency", help = "path to folder with frequency files ", default  = 'data/frequencies/')
parser.add_option("-s", "--start", dest = "start", help = "start position of the region", default = "0")


(options, args) = parser.parse_args()

# Checking if the path ends with '/'
if (not options.ancestor.endswith('/')):
	options.ancestor = options.ancestor+'/'

if (not options.genome.endswith('/')):
	options.genome = options.genome+'/'

if (not options.frequency.endswith('/')):
	options.frequency = options.frequency+'/'

# Checks whenever a input was given for OptionParser.
if len(options.chromosome) == 0:
	sys.exit('No chromosome(s) specified, Program closed prematurely.')
if len(options.start) == 0:
	sys.exit('No start position specified, Program closed prematurely.')


# Create list for ancestor input
anc_list = []
for fn in os.listdir(options.ancestor):
	if fn.endswith('fa'): 
		anc_list.append(fn)
anc_list = natsort.natsorted(anc_list)

# Create list for ancestor input
ref_list = []
for fn in os.listdir(options.genome):
	if fn.endswith('fa'): 
		ref_list.append(fn)
ref_list = natsort.natsorted(ref_list)

# Create list for ancestor input
freq_hi_list = []
freq_lo_list = []
for fn in os.listdir(options.frequency):
	if fn.endswith('out') and "reversed" not in fn: 
		freq_hi_list.append(fn)
	elif fn.endswith('out') and "reversed" in fn:
		freq_lo_list.append(fn)

freq_hi_list = natsort.natsorted(freq_hi_list)
freq_lo_list = natsort.natsorted(freq_lo_list)


# Creates a list for the chromosomes
chr_list = options.chromosome.split(',')

# Check if the nt are written in upper (unmasked) or lower case (masked), 
# depending on this, the script writes it in one of the two output files.
# change this later as well maybe or move it down
def write_line(chrom, pos, reference, an, output_low, output_high):
	if (reference in 'ACGTacgt') and (an in 'ACGTacgt'):
		if an.islower():
			output_low.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n" %(chrom, pos + 1, reference, an))
		else:
			output_high.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n" %(chrom, pos + 1, reference, an))


for chromosome in chr_list:
	print(chromosome)

	ancestor_file = [x for x in anc_list if chromosome in x][0]
	print(ancestor_file)
	reference_file = [x for x in ref_list if chromosome in x][0]
	print(reference_file)
	frequency_file = [x for x in freq_hi_list if chromosome in x][0]
	print(frequency_file)

	# Define input files for the ancestral sequence and the reference. 
	ancestor_fasta = pysam.Fastafile(options.ancestor + ancestor_file)
	ref_fasta = pysam.Fastafile(options.genome + reference_file)

	# Open population frequency files
	# Note: frequency files are not vcf files
	input_snps = open(options.frequency + frequency_file, 'r')

	# Check if the ancestor seq is of same length as the reference
	if ancestor_fasta.nreferences != 1:
		sys.exit('There are more than 1 record in the fasta file.')
	else:
		ancestor_record = ancestor_fasta.references[0]
		if chromosome in ancestor_record:
			anc_length = ancestor_fasta.get_reference_length(ancestor_fasta.references[0])
		elif 'chr' + chromosome in ancestor_record:
			anc_length = ancestor_fasta.get_reference_length(ancestor_fasta.references[0])
		else:
			sys.exit('The requested chromosome cannot be found in the record of the ancestor fasta\n%s\n%s'%('%s:'%chromosome, ancestor_record))

	if anc_length != ref_fasta.get_reference_length(ref_fasta.references[0]):
		sys.exit('Ancestor fasta and ref fasta have not the same length, Chromosome: %s'%('%s'%chromosome))
	else:
		print('Ancestor sequence control: chr ' + str(chromosome) +', sequence is good')

	# Create output files for upper and lower case variants.
	output_high = open("derived_variants_chr%s_case_upper.vcf" %chromosome, 'w')
	output_low = open("derived_variants_chr%s_case_lower.vcf" %chromosome, 'w')
	output_low.write("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	output_high.write("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

	# Split the ancestor and reference sequences in chunks of 350000.
	# Done so that the script reads this in chunks. 
	index_outer = (int(options.start), int(options.start) + 350000)
	snp_line = input_snps.readline()
	snp_line = input_snps.readline()
	snp_pos = int(snp_line.strip().split('\t')[1])

	while index_outer[0] <= anc_length:
		# Fetch the reference and ancestral sequence.
		ref_seq = ref_fasta.fetch(ref_fasta.references[0], index_outer[0], index_outer[1])
		anc_seq = ancestor_fasta.fetch(ancestor_record, index_outer[0], index_outer[1])
			
		# Itterate over the ancestral sequence position. 
		for i, anc in enumerate(anc_seq):
				
			# Enumerate sets i always to 0, therefore it needs to be summed with the starting position of index_outer
			# Checking if current position in the sequences is equal to the snp position in the frequency files
			if i + index_outer[0] == (snp_pos-1):
				snp_allele = snp_line.strip().split('\t')[5].split(':')[0]
					
				# Checking if Case 2 is true, i.e. ancestor and reference are same, but there is a different allele in pop at >90%. 
				if (anc == ref_seq[i]) and (anc != snp_allele) and (anc != '-') and (anc != '.'):
					# Case 2 write out
					write_line(chromosome, i + index_outer[0], ref_seq[i], snp_allele, output_low, output_high)
					# 'Try' is used, to avoid error when the end of snp file is reached but not yet the end of the sequence. 
					try: 
						snp_line = input_snps.readline()
						snp_pos = int(snp_line.strip().split('\t')[1])
					except:
						pass
 
				# Only take next line in snp file
				# so just continue to next if it does not match anything else?
				else: 
					try:
						snp_line = input_snps.readline()
						snp_pos = int(snp_line.strip().split('\t')[1])
					except: 
						pass

			elif anc != '-' :

				# Checking if Case 1 is true: i.e. ancestor and reference are different. 
				if (anc != ref_seq[i]):
					# Case 1 write out
					write_line(chromosome, i + index_outer[0], ref_seq[i], anc, output_low, output_high)
					# Open population frequency files

		index_outer = (index_outer[1], index_outer[1] + 350000)

	# Go through the list of SNPs in the other file where the alternative is the ancestral
	# and the reference is >90% i.e. where case 3 is true. 
	frequency_file = [x for x in freq_lo_list if chromosome in x][0]
	print(frequency_file)

	# Note: frequency files are not vcf files
	input_snps = open(options.frequency + frequency_file, 'r')

	# Split the ancestor and reference sequences in chunks of 350000.
	# Done so that the script reads this in chunks. 
	index_outer = (int(options.start), int(options.start) + 350000)
	snp_line = input_snps.readline()
	snp_line = input_snps.readline()
	snp_pos = int(snp_line.strip().split('\t')[1])

	while index_outer[0] <= anc_length:
		# Fetch the reference and ancestral sequence.
		ref_seq = ref_fasta.fetch(ref_fasta.references[0], index_outer[0], index_outer[1])
		anc_seq = ancestor_fasta.fetch(ancestor_record, index_outer[0], index_outer[1])

		# Checking if case 3 is true
		if (anc != ref_seq[i]) and (anc == snp_allele) and (anc != '-') and (anc != '.'):
			# Case 3 write out
			write_line(chromosome, i + index_outer[0], ref_seq[i], snp_allele, output_low, output_high)
			try: 
				snp_line = input_snps.readline()
				snp_pos = int(snp_line.strip().split('\t')[1])
			except: 
				pass
				
		index_outer = (index_outer[1], index_outer[1] + 350000)		
				
output_high.close()
output_low.close()


