#!/usr/bin/env python
'''
:Author: Seyan Hu
:Date: 10-11-2022
:Usage: python <script.py> -s <path to simulated variants> -d <path to derived variants>

Trims the vcf file of simulated variants to the same amount of derived variants. 

:Example:
python trim_simulated_var.py -s ./ -d ./
'''

# Import dependencies.
import sys, os
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-s", "--simu", dest="simu", help="Simulated variants", default= './') 
parser.add_option("-d", "--derived", dest="derived", help="Derived variants", default= './') 

(options, args) = parser.parse_args()


## FIX THIS; DO WE NEED TO CHOP SIMULATED OR CONCATENATE DERIVED?
# Sorted list of files in directory
simu_file_list = []
derived_file_list = []
for fn in os.listdir(options.simu):
	if fn.startswith('snps_sim'): # simu
		simu_file_list.append(fn)
#for fn in os.listdir(options.derived):
#	if fn.startswith('snps_derived_') and fn.endswith('_upper.vcf'): # derived
#		derived_file_list.append(fn)
		
simu_file_list = sorted(simu_file_list)
#derived_file_list = sorted(derived_file_list)

#print(simu_file_list)
#print(derived_file_list)


# Iterate over list of files. 
for simu_fn in simu_file_list:
	#derived_fn = derived_file_list[num]
	
	### THIS DOES NOT WORK ANYLONGER
	### HARDCODED NOW
	num = simu_fn.split('_site_')[-1]
	num = num.split('.')[0]
	derived_fn = 'snps_derived_var_chr_' + str(1) + '_case_upper.vcf'
	
	#print(simu_fn, derived_fn)

	# Iterate over lines in derived and simulated variant files.
	# Count the lines in derived vcf file. (2 header lines)
	d_lines_count = 0
	for d_line in open(options.derived + derived_fn, 'r'):
		if not '#' in d_line:
			d_lines_count += 1


	# Write simulated lines to output, 
	# Removes headers and RANDOM sample varaints to new file (should have same number of variants as derived).
	os.system("sed -i.bak '/#/d' " + options.simu + simu_fn)
	os.system('gshuf -n '+str(d_lines_count)+' '+options.simu + simu_fn+' > '+'trimmed_temp_' + simu_fn)
	## sh: shuf: command not found
	## brew install coreutils


	# Add headers to file and sort on position. 
	#os.system('echo "##fileformat=VCFv4.1" > '+'trimmed_t2_' + simu_fn)
	os.system('echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" > '+'trimmed_t2_' + simu_fn)
	os.system('cat trimmed_temp_' + simu_fn+' >> '+'trimmed_t2_' + simu_fn)
	os.system('rm trimmed_temp_' + simu_fn)
	os.system("sort -t'\t' -k2g " + 'trimmed_t2_' + simu_fn + " -o " + 'trimmed_t3_' + simu_fn)
	
	os.system('echo "##fileformat=VCFv4.1" > '+'trimmed_' + simu_fn)
	os.system('cat trimmed_t3_' + simu_fn+' >> '+'trimmed_' + simu_fn)
	
	os.system('rm trimmed_t2_' + simu_fn)
	os.system('rm trimmed_t3_' + simu_fn)
	os.system('rm '+options.simu + simu_fn+'.bak')
	
	


# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_simulated_vcf_trimming.txt', 'x')
indication.close
os.rename('./finished_simulated_vcf_trimming.txt', './output/finished_simulated_vcf_trimming.txt')






