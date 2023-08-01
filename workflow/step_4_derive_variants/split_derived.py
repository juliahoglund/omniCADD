"""
:Author: Christian Gross
:Contact: cgross@tudelft.nl
:Date: 01-08-2018

This script takes the vcf file with derived mouse single variations 
and removes variants which are in a sequence of positions and therefore might be longer variants than just SNPs.


:Edited by: Seyan Hu
:Date: 20-10-2022
:Usage: python <python file> <path to one vcf files> <file identifier>

:Example:
python filtering_derived_for_singletons.py ./ derived_var_chr_19_case_

"""


# Import dependencies
import os,sys
from optparse import OptionParser


# change this as well to something better
# Assign input to variables
path_vcf = sys.argv[1]		# './'
file_ident = sys.argv[2]	# 'derived_var_chr_1_case_'


# OptionParser for the previously generated derived variants files (vcf)
parser = OptionParser()
#parser.add_option("-i", "--input", dest="input", help="Path to derived variants.", default='./derived_var_chr_19_case_lower.vcf') 
parser.add_option("-p", "--path", dest="path", help="Path to folder with vcf files.", default= path_vcf) 
(options, args) = parser.parse_args()


# Error query if the path is empty  
if options.path == '':
	sys.exist('Path to model files is not given.')


# Checking if the path ends with '/' 
if (not options.path.endswith('/')):
	options.path = options.path+'/'


# Makes a list from the files in the directory and remove 'other' from files
file_list = []
for file in os.listdir(options.path):
	if file.startswith(file_ident):
		file_list.append(file)
#print(file_list)


## CHECK SO IT ALL MAKES SENSE; ALSO DOES IT NEED A WRAPPER?
# Iterate over files in list. 
for file in file_list:
	
	
	# Read vcf file and create two seperate files for snps and indels.
	infile = open(options.path + file,'r')
	outfile1 = open(options.path + 'snps_' + file,'w')
	outfile2 = open(options.path + 'indels_' + file,'w')

	outfile1.write('##fileformat=VCFv4.1\n')
	outfile1.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')

	outfile2.write('##fileformat=VCFv4.1\n')
	outfile2.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')


	# Loop through the lines in vcf file
	counts = 0
	previous_pos = 0
	previous_lines = ''
	jump = False
	for lines in infile:
		
		# Skips header
		if lines.startswith('#'):
			continue
		
		# Looks through lines containing information
		else:
			line = lines.split('\t')
			# If 'previous_line' is an empty string, it appends the current line to the previous line 
			# and sets 'previous_pos' to the variants position of the current line.
			# Happens when the loop encounters the first line of information. 
			if previous_lines == '':
				previous_lines = lines
				previous_pos = int(line[1])
				continue
			# If the position of the current variant is next to the previous position, 
			# it writes the line to the indel file
			# and sets the 'previous_line' to the current line 
			# and 'previous_pos' to the variants position of the current line.
			# It also sets 'jump' to True.
			if (int(line[1])-1==previous_pos) or (int(line[1])==previous_pos):
				jump = True
				outfile2.write(previous_lines)
				previous_lines = lines
				previous_pos = int(line[1])
				continue
			# If the position of the current variant is not next to the previous position and not the first line,
			# it check whenever 'jump' is still True. 
			# If 'jump' is still True, it will be set to False and the line is written to the indel file. 
			# The 'previous_line' is set to the current line 
			# and 'previous_pos' to the variants position of the current line.
			
			# Else it writes the line to the snp file and 
			# the 'previous_line' is set to the current line 
			# and 'previous_pos' to the variants position of the current line.
			else:
				if jump == True:
					jump = False
					outfile2.write(previous_lines)
					previous_lines = lines
					previous_pos = int(line[1])
					continue
				else:
					outfile1.write(previous_lines)
					previous_lines = lines
					previous_pos = int(line[1])

	# If 'jump' is still equal to True write the lines to the indel file,
	# else write it to the snp file. 
	if jump == True:
		outfile2.write(previous_lines)
	else:
		outfile1.write(previous_lines)


# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_snp_filtering.txt', 'x')
indication.close()
os.rename('./finished_snp_filtering.txt', './output/finished_snp_filtering.txt')


