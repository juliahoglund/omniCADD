#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 31-10-2022
:Extension and modification: Julia Höglund
:Date: 02-05-2023

:Usage: python <script.py> -a <ancestor seq file> -r <reference seq file> -c <chromosome number> -o <output file name>

:Example:
python create_parameters.py -a Ancestor_1.fa -r chr1.fa -c 1 -o chr1.log

Gathers parameters for the generation of the simulated variants for one chromosome.
This script is performed for one chromosome thus a wrapper is needed.
It gathers the counts of nucleotide substitutions, insertions and deletions.
This script takes in one chromosome, and for severeal, has to be called with the wrapper,
wrapper_create_parameters.py

For this script to work properly, the input fasta files had to be single stranded.


:Sanity check:
anc = 'GAAGCCGTGGAGCAGGG----------------------------AGCACCAGCGGGCGACGGTGGAGGACATGGG--------GCTGGCCGGGCAGAGG-ACGCAAA'
ref = 'AAAGCTGCAGATCGGACGGGTGAGAGGGTGGGGGTCAGAGGTCAGAGCGTTTACGGGCGACGGTGGAGGACATAGGCCGGGGTAGCTGGTGGGGCTGCGANNCGANNA'

22 totalrefA 	108 total nt count		36 nt insertions		1 ACn		1 TCn		1 CA		0 CAn		1 GA		5 GAn
17 totalrefC 	27 total CpG count		3 nt deletions			2 AGn		0 TGn		1 CG		0 CGn		0 GC		1 GCn
50 totalrefG 	14 mutations			1 gap in both			2 ATn		0 TAn		2 CT		2 CTn		1 GT		0 GTn
15 totalrefT	6 CpG mutations

Insertion: 	length = 28, count = 1;		length = 8, count = 1
Deletion:	length = 1, count = 1; 		lenght = 2, count = 1

"""

# Import dependencies
import os, sys
from optparse import OptionParser
from collections import defaultdict
from progressbar import Percentage, ProgressBar, Bar, ETA # allow progress bar

# OptionParser for input.
parser = OptionParser()
parser.add_option("-a", "--ancestor", dest="ancestor", help="extracted ancestor sequence",default="output/extracted_ancestor/Ancestor_Pig_Cow.1_chr1.fa")
parser.add_option("-r", "--reference", dest="reference", help="reference sequence of the species of interest",default="genome/")
parser.add_option("-c", "--chr", dest="chr", help="chromosome from which variants should be simulated",default="1")
parser.add_option("-o", "--outfile", dest="outfile", help="name of output file",default= 'chr1.log')

(options, args) = parser.parse_args()


# Open inputs.
anc_open = open(options.ancestor)
anc_open.readline() # empty
anc_open.readline() # header
anc_str = anc_open.readline().replace('\n', '') # sequence

ref_open = open(options.reference)
ref_open.readline() # empty
ref_open.readline() # header
ref_str = ref_open.readline().replace('\n', '') # sequence

# Function for turning seq into string and then to list per nt.
def stringTOlist(seq):
	list_nt = [nt for nt in seq]
	return list_nt

list_anc = stringTOlist(anc_str)
list_ref = stringTOlist(ref_str)

# Total counts and mutation counts.
total = 0
mut = 0
totalCpG = 0
mutCpG = 0
# Total nt in reference.
totalrefA = 0
totalrefC = 0
totalrefG = 0
totalrefT = 0
# Total nt substitution counts.
ACn = 0
AGn = 0
ATn = 0
CAn = 0
CGn = 0
CTn = 0
GAn = 0
GCn = 0
GTn = 0
TAn = 0
TCn = 0
TGn = 0
# Total nt substitution counts in CpG sites.
CA = 0
CG = 0
CT = 0
GA = 0
GC = 0
GT = 0
# Insertions and deletions counts.
nt_insertion = 0
nt_deletion = 0
gap = 0
dict_insertion = {}
dict_deletion = {}


# Calc total nt.
len_seq = len(list_anc)
total += len_seq


print('Start counting! Chr: ' + options.chr)

# Iterate through nt and compare.
N = len(list_ref)
pbar = ProgressBar(widgets=[Bar('=', '[', ']'), ' ', Percentage(), ' ', ETA()], maxval = N).start()

for num, nt_ref in enumerate(list_ref):
	nt_anc = list_anc[num]

	## Checks the nt from the reference counts.
	if nt_ref == 'G':
		totalrefG += 1
	elif nt_ref == 'C':
		totalrefC += 1
	elif nt_ref == 'A':
		totalrefA += 1
	elif nt_ref == 'T':
		totalrefT += 1

	## Check for insertion, gaps, deletion counts.
	# Append the 2 nt (ancestor-reference).
	list_nt_anc_ref = nt_anc + nt_ref
	# Check for insertions (when the nt is '-' (a gap) in ancestor but in ref a C, G, A or T).
	if list_nt_anc_ref.startswith('-') and not list_nt_anc_ref.endswith('N'):
		nt_insertion += 1
	# Check for gaps in both sequences.
	elif list_nt_anc_ref.startswith('-') and list_nt_anc_ref.endswith('N'):
		gap += 1
	# Checks for deletions (when the nt is not '-' (so a C, G, A or T) in ancestor but in ref a N).
	elif not list_nt_anc_ref.startswith('-') and list_nt_anc_ref.endswith('N'):
		nt_deletion +=1

	## Check for substitution counts.
	# Checks if ancestor nt starts with 'A'.
	elif list_nt_anc_ref.startswith('A'):
		if list_nt_anc_ref.endswith('C'):
			ACn += 1
			mut += 1
		elif list_nt_anc_ref.endswith('G'):
			AGn += 1
			mut += 1
		elif list_nt_anc_ref.endswith('T'):
			ATn += 1
			mut += 1
	# Checks if ancestor nt starts with 'T'.
	elif list_nt_anc_ref.startswith('T'):
		if list_nt_anc_ref.endswith('C'):
			TCn += 1
			mut += 1
		elif list_nt_anc_ref.endswith('G'):
			TGn += 1
			mut += 1
		elif list_nt_anc_ref.endswith('A'):
			TAn += 1
			mut += 1

	## Check for CpG/GpC sites folowed by counting and for substitution counts for C and G to another nt.
	# Check if ancestor nt starts with 'C'.
	elif list_nt_anc_ref.startswith('C'):
		# Checks if it is a CpG/GpC site.
		# Does not count as normal substitution.
		if (num > 0 and list_anc[num - 1] == 'G') or (num < (len_seq - 1) and list_anc[num + 1] == 'G'):
			# Counts CpG/GpC sites occurences.
			totalCpG += 1
			# Checks for mutations in the site.
			if list_nt_anc_ref.endswith('A'):
				CA += 1
				mutCpG += 1
			elif list_nt_anc_ref.endswith('G'):
				CG += 1
				mutCpG += 1
			elif list_nt_anc_ref.endswith('T'):
				CT += 1
				mutCpG += 1
		# If not in site, mutations are still checked.
		else:
			if list_nt_anc_ref.endswith('A'):
				CAn += 1
				mut += 1
			elif list_nt_anc_ref.endswith('G'):
				CGn += 1
				mut += 1
			elif list_nt_anc_ref.endswith('T'):
				CTn += 1
				mut += 1
	# Check if ancestor nt starts with 'G'.
	elif list_nt_anc_ref.startswith('G'):
		# Checks if it is a CpG/GpC site.
		if (num > 0 and list_anc[num - 1] == 'C') or (num < (len_seq - 1) and list_anc[num + 1] == 'C'):
			# Counts CpG/GpC sites occurences.				totalCpG += 1
			# Checks for mutations in the site.
			if list_nt_anc_ref.endswith('A'):
				GA += 1
				mutCpG += 1
			elif list_nt_anc_ref.endswith('C'):
				GC += 1
				mutCpG += 1
			elif list_nt_anc_ref.endswith('T'):
				GT += 1
				mutCpG += 1
		# If not in site, mutations are still checked.
		else:
			if list_nt_anc_ref.endswith('A'):
				GAn += 1
				mut += 1
			elif list_nt_anc_ref.endswith('C'):
				GCn += 1
				mut += 1
			elif list_nt_anc_ref.endswith('T'):
				GTn += 1
				mut += 1
	pbar.update(num)
pbar.finish()

## Checks for insertions and deletions.
# Dictionaries for instertions and deletions (length(key), count(value)).
insertionsizes = defaultdict(int)
insertion = 0
deletionsizes = defaultdict(int)
deletion = 0

print("Measuring indels ...")
# Iterate over the nt in the sequence.
N = len_seq
pbar = ProgressBar(widgets=[Bar('=', '[', ']'), ' ', Percentage(), ' ', ETA()], maxval = N).start()

for position in range(1, len_seq + 1):
	nt_anc = list_anc[position-1]
	nt_ref = list_ref[position-1]

	# If there is an nt insertion found add 1 to insertion (size), if it is the start of insertion or previous nt was also an insertion.
	if nt_anc == "-" and ((insertion > 0) or (insertion == 0 and nt_ref in ["A","C","G","T"])):
		insertion += 1
	# Else (at end of insertion) add the insertion size to the dict and add a count to it.
	# Afterwards set the size ('insertion') back to 0.
	elif nt_ref in ["A","C","G","T"]:
		if insertion > 0:
			insertionsizes[insertion] += 1
			insertion = 0

	# If there is an nt deletion found add 1 to deletion (size), if it is the start of deletion or previous nt was also a deletion.
	if nt_ref == "N" and ((deletion > 0) or (deletion == 0 and nt_anc in ["A","C","G","T"])):
		deletion += 1
	# At the end of the deletion add the size to the dict and add a count to it.
	# Afterwards set the size back to 0.
	elif nt_anc in ["A","C","G","T"]:
		if deletion > 0:
			deletionsizes[deletion] += 1
			deletion = 0
	pbar.update(position)
pbar.finish()

# Define controls.
chrom, window_start, pos = options.chr, 1, position
cmut = mut
ctotal = total
cmutCpG = mutCpG
ctotalCpG = totalCpG
cA, cC, cG, cT = totalrefA, totalrefC, totalrefG, totalrefT


anc_open.close()
ref_open.close()
# Create output.
'''
:Output format:
###CHROM <num>
chrom	window_start	pos	cmut	ctotal	cmutCpG	ctotalCpG	cA	cC	cG	cT
##STATS
#A	C	G	T	CpGs
...
#y	N	AC	AG	AT	CA	CG	CT	GA	GC	GT	TA	TC	TG
...
#yCpG	NCpG	CA	CG	CT	GA	GC	GT
...
##INSERTIONS
#len	count
...
##DELETIONS
#len	count
...
'''

print("Writing output.")
output = open(options.outfile,'w')

# Write to output.
output.write("###CHROM " + str(options.chr) + "\n")
output.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"%(chrom, window_start, pos, cmut, ctotal, cmutCpG, ctotalCpG, cA, cC, cG, cT))

output.write("\n##STATS\n")
output.write("#A\tC\tG\tT\tCpGs\n")
output.write("\t".join(map(str,[totalrefA, totalrefC, totalrefG, totalrefT, totalCpG])) + "\n")
output.write("#y\tN\tAC\tAG\tAT\tCA\tCG\tCT\tGA\tGC\tGT\tTA\tTC\tTG\n")
output.write("\t".join(map(str,[mut, total, ACn, AGn, ATn, CAn, CGn, CTn, GAn, GCn, GTn, TAn, TCn, TGn])) + "\n")
output.write("#yCpG\tNCpG\tCA\tCG\tCT\tGA\tGC\tGT\n")
output.write("\t".join(map(str,[mutCpG, totalCpG, CA, CG, CT, GA, GC, GT])) + "\n")

output.write("##INSERTIONS\n")
output.write("#len\tcount\n")
for key, value in insertionsizes.items():
	output.write("%d\t%d\n"%(key, value))

output.write("##DELETIONS\n")
output.write("#len\tcount\n")
for key, value in deletionsizes.items():
	output.write("%d\t%d\n"%(key, value))

print("Done!")
# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_create_parameters.txt', 'x')
indication.close()
os.rename('./finished_create_parameters.txt', './output/finished_create_parameters.txt')
