#!/usr/bin/env python
# -*- coding: ASCII -*-

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
:Usage: python derive_variants.py.py -c <chr num> -a <path to ancestor seq> -g <path to genome> -f <path to frequency files> -s <start position of region>

:Example:
python derive_variants.py -c 1,2,3,4,5 -a data/extracted_ancestor -g data/genome -f data/frequencies -s 0
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

freq_list = []
# Create list for ancestor input
for fn in os.listdir(options.frequency):
	if fn.endswith('out'): 
		freq_list.append(fn)

freq_list = natsort.natsorted(freq_list)

# Creates a list for the chromosomes
chr_list = options.chromosome.split(',')

# Check if the nt are written in upper (unmasked) or lower case (masked), 
# depending on this, the script writes it in one of the two output files.
# change this later as well maybe or move it down
def write_line(chrom, pos, reference, an, output):
	if (reference in 'ACGTacgt') and (an in 'ACGTacgt'):
		if an.islower():
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n" %(chrom, pos + 1, reference, an.upper()))
		else:
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n" %(chrom, pos + 1, reference, an))


for chromosome in chr_list:
	print(chromosome)

	ancestor_file = [x for x in anc_list if chromosome in x][0]
	print(ancestor_file)
	reference_file = [x for x in ref_list if chromosome in x][0]
	print(reference_file)
	frequency_file = [x for x in freq_list if chromosome in x][0]
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
	output = open("derived_variants_chr%s.vcf" %chromosome, 'w')
	output.write("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

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
			
			# Enumerate sets i always to 0, therefore i needs to be summed with the starting position of index_outer
			# Checking if current position in the sequences is equal to the snp position in the frequency files
			if i + index_outer[0] == (snp_pos-1):
				snp_allele = snp_line.strip().split('\t')[5].split(':')[0]
				
				# Checking if Case 1 derived snp is given
				if (anc == ref_seq[i]) and (anc != snp_allele) and (anc != '-') and (anc != '.'):
					# Case 1 write out
					print("case1")
					write_line(options.chrom, i + index_outer[0], ref_seq[i], snp_allele, output_low, output_high)
					#print i+index_outer[0],anc,ref_seq[i],snp_allele
					# 'Try' is used, to avoid error when the end of snp file is reached but not yet the end of the sequence. 
					try: 
						snp_line = input_snps.readline()
						snp_pos = int(snp_line.strip().split('\t')[1])
					except:
						pass
				
				# Checking if case 3 (or 2) derived snp is given
				elif (anc != ref_seq[i]) and (anc != snp_allele) and (anc != '-'):
					print("case2,3")
					# Case 3 write out
					write_line(options.chrom, i + index_outer[0], ref_seq[i], snp_allele, output_low, output_high)
					#print i+index_outer[0],ref_seq[i],snp_allele
					try: 
						snp_line = input_snps.readline()
						snp_pos = int(snp_line.strip().split('\t')[1])
					except: 
						pass
				
				# Only take next line in snp file
				else: 
					print("here") 
					try:
						snp_line = input_snps.readline()
						snp_pos = int(snp_line.strip().split('\t')[1])
					except: 
						pass
						
			elif anc != '-' :
				print("anc is ref if not followed by ref is not")
				# Checking if case 5 (or 3) derived variant
				if (anc != ref_seq[i]):
					print("ref is not anc:"+ref_seq[i]+anc)
					# Case 5 write out
					#print i+index_outer[0],anc
					write_line(options.chrom, i + index_outer[0], ref_seq[i], anc, output_low, output_high)
					
		index_outer = (index_outer[1], index_outer[1] + 350000)		
				
output.close()


