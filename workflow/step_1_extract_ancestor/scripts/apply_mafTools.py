#!/usr/bin/env python

"""
:Author: Seyan Hu
:Date: 28-9-2022
:Extension and modification: Julia Höglund
:Date 2023-03-17

:Usage: python apply_mafTools.py -m <path to marked files> -g <path to reference genome> -f <path to filtered files> -o [ref_species,ancestor,out_species] -s <path to flipped files>

:mafDuplicateFilter:
Removes all duplicated sequences and keeps only
the one sequence that is the most similar to the block consensus
:mafStrander:
Flips all alignment blocks in which the species of intrest and its
ancestors have been on the negative strand
:mafRowOrderer:
Reorders species within any alignment block

:Example:
python apply_mafTools.py -m ./data/marked -g ./genome/chr -f mDF/mDF_* -o Sus_scrofa,Ancestor_Pig_Cow_Sheep,Bos_tauros -s mStr/mStr_

"""

# Import dependencies
import sys, os, subprocess, argparse
import gzip
from sh import gunzip
from optparse import OptionParser, OptionGroup

# OptionParser for the input directory (process duplicates output), to be marked ancestor seq, and file identifier.
parser = OptionParser()
group = OptionGroup(parser, "Usage (extended)", "This script takes as input the maf files where the \
            ancestor of choice has been marked in the previous step. \
			The input is given like this: \
						python apply_mafTools.py -m <path to marked files> -g <path to reference genome> -f <path to filtered files> -o [ref_species,ancestor,out_species] -s <path to flipped files> \
    Where: \
                                                            <PATH TO MARKED FILES> is the path to the ancestor marked msa files (maf, gzipped or standard) to be further analysed \
                 <PATH TO REFERENCE GENOME> is the path to the folder with the reference genome of choice, to check stranding \
                    <PATH TO FILTERED FILES> is the path (including file prefix) that will be created during the first step, mafDuplicateFilter \
[REF_SPECIES,ANCESTOR,OUT_SPECIES] is the names of the species (separated by a comma, no space) to be ordered, reference species, the marked ancestor and the out species all that were chosen and produced with mark_ancestor.py \
                                            <PATH TO FLIPPED FILES>is the path (including file prefix) that will be created during the second step, mafStrander")


parser.add_option_group(group)
# mafDuplicateFilter
parser.add_option("-m", "--marked", dest="marked", help="path to folder (including prefix of files) with marked alignment files", default= "./marked/marked_")
# mafStrander
parser.add_option("-g", "--genome", dest="genome", help="path to (reference) genome including the prefix", default="./")
parser.add_option("-f", "--filter", dest="filtered", help="path to folder (including prefix of files) with filtered files (to be computed)", default= "mDF/mDF_")
# mafRowOrderer
parser.add_option("-o", "--order", dest="order", help="species and ancestor in the desired order, separated only by commas", default="")
parser.add_option("-s", "--stranded", dest="stranded", help="path to folder (including prefic of files) with files with flipped strands", default="mStr/mStr_")
parser.add_option("-r", "--rowpath", dest="rowpath", help="prefix to use for the last files, after reordering rows", default="mRO/mRO_")
# clean intermediate files
parser.add_option("-c", "--clean", dest="clean", help="remove intermediate files after computation [yes/no]", default="no")
parser.add_option("-p", "--previous", dest="previous", help="remove marked files from previous step (mark ancestor) after computation [yes/no]", default="no")

if len(sys.argv)==1:
    parser.print_help()
    parser.exit()
(options, args) = parser.parse_args()

# CHECK input
if options.marked == '':
    sys.stderr.write("Found no files with given path and identifier (prefix).\n")
    sys.exit()
if options.genome == '':
    sys.stderr.write("Found no path matching the input.\n")
    sys.exit()
if options.order == '':
    sys.stderr.write("No names / erroneously named species given!\n")
    sys.exit()

##############################
##### mafDuplicateFilter #####
##############################
file_prefix = options.marked.split('/')[-1]
path = '/'.join(options.marked.split('/')[0:-1]) + '/'

filenames = [x for x in os.listdir(path) if x.startswith(file_prefix)]

filter_file_prefix = options.filtered.split('/')[-1]
filter_path = '/'.join(options.filtered.split('/')[0:-1]) + '/'

## if path not exist makedir mDF else not do it
if not os.path.isdir(filter_path):
    os.system('mkdir ' + str(filter_path))

# Loop through files
for file in filenames:

    ## open file if gzipped
    if file.endswith(".gz"):
        gunzip(str(path)+str(file))

    file = file.replace('.gz', '')
    print('Performing DuplicateFilter for: ' + file)
	# Perform mafTools on commandline
    os.system('mafDuplicateFilter --maf ' + str(path) + str(file) + ' > ' + str(filter_path) + str(filter_file_prefix) + file)
    os.system('gzip ' + str(path) + str(file))
    os.system('gzip ' + str(filter_path) + str(filter_file_prefix) + file)

#######################
##### mafStrander #####
#######################

file_prefix = options.filtered.split('/')[-1]
path = "/".join(options.filtered.split('/')[0:-1]) + '/'

filenames = [x for x in os.listdir(path) if x.startswith(file_prefix)]

stranded_file_prefix = options.stranded.split('/')[-1]
stranded_path = '/'.join(options.stranded.split('/')[0:-1]) + '/'

## if path not exist makedir mDF else not do it
if not os.path.isdir(stranded_path):
    os.system('mkdir ' + str(stranded_path))

# Loop through files
for file in filenames:

    ## open file if gzipped
    if file.endswith(".gz"):
        gunzip(str(path)+str(file))

    file = file.replace('.gz', '')
    # Split path from the file name itself.
    file = file.split('/')[-1]
    print('Performing Strander for: ' + file)

	# Get chr number.
    chr_num = file.split('epo.')[-1] # in case all is not epo change?
    chr_num = chr_num.split('_')[0]
    seq_inp = options.genome + '.' + chr_num

	# Perform mafTools on commandline
    os.system('mafStrander --maf ' + str(path) + str(file) + ' --seq ' + seq_inp + ' --strand + > ' + str(stranded_path) + str(stranded_file_prefix) + file)

    if options.clean=='yes':
        os.system('rm ' + str(path) + str(file))
    else:
        os.system('gzip ' + str(path) + str(file))
    os.system('gzip ' + str(stranded_path) + str(stranded_file_prefix) + file)

#########################
##### mafRowOrderer #####
#########################

file_prefix = options.stranded.split('/')[-1]
path = "/".join(options.stranded.split('/')[0:-1]) + '/'

ordered_file_prefix = options.rowpath.split('/')[-1]
ordered_path = '/'.join(options.rowpath.split('/')[0:-1]) + '/'

## if path not exist makedir mDF else not do it
if not os.path.isdir(ordered_path):
    os.system('mkdir ' + str(ordered_path))

filenames = [x for x in os.listdir(path) if x.startswith(file_prefix)]

# Loop through files
for file in filenames:

    ## open file if gzipped
    if file.endswith(".gz"):
        gunzip(str(path)+str(file))

    file = file.replace('.gz', '')
    print('Performing RowOrderer for: ' + file)
	# Perform mafTools on commandline
    os.system('mafRowOrderer --maf ' + str(path) + str(file) + ' --order ' + str(options.order) + ' > ' + str(options.rowpath) + file)

    if options.clean=='yes':
        os.system('rm ' + str(path) + str(file))
    else:
        os.system('gzip ' + str(path) + str(file))
    os.system('gzip ' + str(options.rowpath) + str(file))

if options.clean=='yes':
    os.system('rm -r ' + stranded_path + ' ' + filter_path)

if options.previous=='yes':
    os.system('rm -r ' + '/'.join(options.marked.split('/')[:-1]))

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_apply_mafTools.txt', 'x')
indication.close
os.rename('./finished_apply_mafTools.txt', './output/finished_apply_mafTools.txt')
