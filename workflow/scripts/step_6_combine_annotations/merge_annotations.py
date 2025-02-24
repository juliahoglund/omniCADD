#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
script to merge annotations to one annotation file per chromosome.
"""

import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v", "--vep",
    help="Processed file (chromosome wide) with vep annotations", 
    type=str, 
    required=True)
parser.add_argument("-b", "--bed",
    help="Bed file (chromosome wide) with processed phast and gerp annotations",
    type=str, 
    required=True)
parser.add_argument("-o", "--outfile",
    help="Name of outfile with combined annotations",
    type=str, 
    required=True)

args = parser.parse_args()

 # combine vep and evolutionary constraint
vepfile = pd.read_csv(args.vep, 
 	sep = "\t",
 	low_memory = False)

bedfile = pd.read_csv(args.bed,
 	sep = " ",
    low_memory = False)

# 1. remove some unwanted columns in bed file, 
#	like chrom end maybe more
bedfile = bedfile.drop(columns = ['chr', 'end'])
bedfile = bedfile.rename(columns={"start": "Pos"})
# 2. left outer join with vep versus bed
left_merged = pd.merge(vepfile, bedfile, how="left", on=["Pos"])
#	save output. 

## add more annotations here later

# 3. write to file
left_merged.to_csv(args.outfile, index = False, sep = "\t") 
