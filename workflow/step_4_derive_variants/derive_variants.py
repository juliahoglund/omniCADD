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


# change all this to options parser
# Assigns variables to input
chromosome_num = sys.argv[1]		# 1
ancestor_seq_path = sys.argv[2]		# 'extracted_ancestor/'
genome_path = sys.argv[3]			# 'step_1/'
freq_path = sys.argv[4]				# 'step_4/'
# make failsaife with last slash


# Check if the nt are written in upper or lower case, 
# depending on this, the script writes it in one of the two output files.
def write_line(chrom, pos, reference, an, output_low, output_high):
	if (reference in 'ACGTacgt') and (an in 'ACGTacgt'):
		if an.islower():
			output_low.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n" %(chrom, pos + 1, reference, an))
		else:
			output_high.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n" %(chrom, pos + 1, reference, an))


# OptionParser for the chr number and start position. 
parser = OptionParser()
parser.add_option("-c", "--chromosome", dest= "chrom", help= "investigated Chromosome", default= chromosome_num)
parser.add_option("-s", "--start", dest= "start", help= "start position of the region", default= "0")
(options, args) = parser.parse_args()


# Checks whenever a input was given for OptionParser.
if len(options.chrom) == 0:
	sys.exit('No Chromosome specified, Program closed prematurely.')
if len(options.start) == 0:
	sys.exit('No Startposition specified, Program closed prematurely.')


## can i do a list files and extract chromosomes all from file list here as before+?? 
# Define input files for the Ancestral sequence and the reference. 
#ancestor_fasta = pysam.Fastafile("../../../generate_ancestral_seq/output/dir_generated_ancestor_seq/Ancestor_Mouse_Rat_%s.fa" %options.chrom)
ancestor_fasta = pysam.Fastafile(ancestor_seq_path + "Ancestor_Pig_Cow." + str(options.chrom) + '_chr' + str(options.chrom) + '.fa') # change to insert full name in script
ref_fasta = pysam.Fastafile(genome_path + 'Sus_scrofa_ref_' + str(options.chrom) + '.fa') # change to insert full name of file or something, like prefix


# Open population frequency files
# Note: frequency files are not vcf files
input_snps = open(freq_path + "%s_freq.out"%options.chrom, 'r')


# Check if the ancestor seq is of same length as the reference
# am i not missing some indents??
if ancestor_fasta.nreferences != 1:
	sys.exit('There are more than 1 record in the fasta file.')
else:
	ancestor_record = ancestor_fasta.references[0]
	#print(ancestor_record)
	#print('.%s'%options.chrom)
	if '.%s'%options.chrom in ancestor_record:
		anc_length = ancestor_fasta.get_reference_length(ancestor_fasta.references[0])
	elif 'chr%s'%options.chrom in ancestor_record:
		anc_length = ancestor_fasta.get_reference_length(ancestor_fasta.references[0])
	else:
		sys.exit('The requested chromosome cannot be found in the record of the ancestor fasta\n%s\n%s'%('.%s:'%options.chrom, ancestor_record))

if anc_length != ref_fasta.get_reference_length(ref_fasta.references[0]):
	sys.exit('Ancestor fasta and ref fasta have not the same length, Chromosome: %s'%('%s'%options.chrom))

print('Ancestor sequence control: Chr' + str(chromosome_num) +', sequence is good')



# Create output files for upper and lower case variants.
output_high = open("derived_var_chr_%s_case_upper.vcf" %options.chrom, 'w')
output_low = open("derived_var_chr_%s_case_lower.vcf" %options.chrom, 'w')
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
		
		# Enumerate sets i always to 0, therefore i needs to be summed with the starting position of index_outer
		# Checking if current position in the sequences is equal to the snp position in the frequency files
		if i + index_outer[0] == (snp_pos-1):
			snp_allele = snp_line.strip().split('\t')[5].split(':')[0]
			
			# Checking if Case 1 derived snp is given
			if (anc == ref_seq[i]) and (anc != snp_allele) and (anc != '-') and (anc != '.'):
				# Case 1 write out
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
				try:
					snp_line = input_snps.readline()
					snp_pos = int(snp_line.strip().split('\t')[1])
				except: 
					pass
					
		elif anc != '-' :

			# Checking if case 5 (or 3) derived variant
			if (anc != ref_seq[i]):
				# Case 5 write out
				#print i+index_outer[0],anc
				write_line(options.chrom, i + index_outer[0], ref_seq[i], anc, output_low, output_high)
				
	index_outer = (index_outer[1], index_outer[1] + 350000)

output_high.close()
output_low.close()


