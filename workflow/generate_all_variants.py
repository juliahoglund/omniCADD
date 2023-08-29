'''
:Author: Seyan Hu
:Date: 10-1-2023
:Usage: python <script.py> -i <Path to genome files>

This script generates all possible varains at every position for the fasta files in the given path. 
And outputs a CSV file for each chromosome. 
Note: Genome files should have this format for their filenames: chr19.fna

'''

# Import dependencies.
import sys,os
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="path to genome files", default= "genome/")
parser.add_option("-p", "--prefix", dest="prefix", help="desired file prefix, eg species name", default= "Sus_scrofa_ref_")


(options, args) = parser.parse_args()


# Sorted list of files in directory.
input_list = []
for fn in os.listdir(options.input):
	if fn.endswith('.fa'):
		input_list.append(fn)
input_list = sorted(input_list)


# Iterate over sorted list.
print('Iterating over input files')
for fn in input_list:
	print('Working on: ' + fn)
	
	# Determine chr number.
	chr_num = fn.split('_')[-1].replace('.fa', '')
	
	# Create output file.
	output = open(str(options.prefix) + str(chr_num) + "_all_variants.vcf", 'w')
	output.write("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	
	# Open chromosome file.
	open_f = open(options.input + fn, 'r')
	
	# Create one single continues sequence.
	seq = ''
	for line in open_f:
		if not line.startswith('>'):
			seq_line = line.replace('\n', '')
			seq += seq_line
	
	# Iterate over the nucleotides in 'seq'.
	pos = 1
	for nt in seq:
		# Write all possible variants to output. 
		if nt == 'A':
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'A', 'T'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'A', 'C'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'A', 'G'))
		elif nt == 'T':
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'T', 'A'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'T', 'C'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'T', 'G'))
		elif nt == 'C':
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'C', 'T'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'C', 'A'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'C', 'G'))
		elif nt == 'G':
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'G', 'T'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'G', 'C'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'G', 'A'))
		#elif nt == 'N':
		#	output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'N', 'T'))
		#	output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'N', 'C'))
		#	output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'N', 'G'))
		#	output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'N', 'A'))
		
		pos += 1

	os.system('gzip ' + str(options.input) + str(fn))
# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_generating_all_variants.txt', 'x')
indication.close()
os.rename('./finished_generating_all_variants.txt', './output/finished_generating_all_variants.txt')
