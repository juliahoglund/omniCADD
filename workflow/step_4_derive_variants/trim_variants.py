#!/usr/bin/env python
'''
:Author: Seyan Hu
:Date: 10-11-2022
:Usage: python <script.py> -s <path to simulated variants> -d <path to derived variants>
## add julia and me and change the example and usage etc etc

Trims the vcf file of simulated variants to the same amount of derived variants. 

:Example:
python trim_simulated_var.py -s ./ -d ./
'''

# Import dependencies.

import sys, os, random
from optparse import OptionParser

# OptionParser for input. 
parser = OptionParser()
parser.add_option("-s", "--simulated", dest = "simulated", help = "path to simulated variants", default = "output/")
parser.add_option("-d", "--derived", dest="derived", help="path to derived variants", default = "output/")
#where do they end up and i only want snps right so filtered takes indels as well?

(options, args) = parser.parse_args()

# Sorted list of files in directory
simu_file_list = []
derived_file_list = []
for fn in os.listdir(options.simulated):
    ## fix with prefix?
    if fn.startswith('snps_simVariants_') and fn.endswith('_filtered.vcf'):
        simu_file_list.append(fn)

## fix with prefix?
for fn in os.listdir(options.derived):
    if fn.startswith('derived_variants_') and fn.endswith('_upper.vcf'):
        derived_file_list.append(fn)

simu_file_list = sorted(simu_file_list)
derived_file_list = sorted(derived_file_list)

# Count the total number of derived variants across all files
total_derived_variants = 0
for derived_fn in derived_file_list:
    with open(options.derived + derived_fn, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                total_derived_variants += 1
print("total amount of derived variants: " + str(total_derived_variants))

# Shuffle the simulated variants and take as many as the total derived variants
for simu_fn in simu_file_list:
    output = open('trimmed_' + simu_fn, 'w')
    with open(options.simulated + simu_fn, 'r') as f:
        all_lines = f.readlines()
        headers = [line for line in all_lines if line.startswith("#")]
        lines = [line for line in all_lines if not line.startswith("#")]
        if len(lines) < total_derived_variants:    
            print("Error: The simulated file {simu_fn} has fewer variants ({len(lines)}) than the total derived variants ({total_derived_variants}). Exiting.")
            sys.exit(1)

        # Print the first 5 lines before shuffling
        print("Before shuffling:")
        for l in lines[:5]:
            print(l.strip())

        random.shuffle(lines)

        # Print the first 5 lines after shuffling
        print("\nAfter shuffling:")
        for l in lines[:5]:
            print(l.strip())
        print("\n")
        for header in headers:
            output.write(header)

        # Now, take the same number of lines as total_derived_variants
        for line in lines[:total_derived_variants]:
            output.write(line)

        output.close()

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_simulated_vcf_trimming.txt', x)
indication.close()
os.rename('./finished_simulated_vcf_trimming.txt', './output/finished_simulated_vcf_trimming.txt')
