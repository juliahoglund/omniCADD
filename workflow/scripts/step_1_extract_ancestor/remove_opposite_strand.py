#!/usr/bin/env python

"""
:Author: Seyan Hu
:Date: 29-9-2022
:Extension and modification: Julia Höglund
:Date: 20-3-2023
:Usage: python remove_opposite_strand.py -p <directory to maf files> -f <file prefix> -r <new path to new files>

Removes lines with sequences that are on the negative strand.

:Example:
python remove_opposite_strand.py -p ./ -f rmSP_mS_ -r forwardStrandOnly

"""

# Import dependencies.
import sys, os
from optparse import OptionParser
import gzip
from sh import gunzip

# OptionParser for input.
parser = OptionParser()
parser.add_option("-p", "--path", dest="path", help="Path to processed maf files",default= './pruned')
parser.add_option("-f", "--fileprefix", dest="prefix", help="Prefix of processed files", default= 'rmSP_mS_')
parser.add_option("-r", "--removed", dest="removed", help="directory within which the new output will be saved", default='forwardStrandOnly')
parser.add_option("-c", "--clean", dest="clean", help="remove intermediate files after computation [yes/no]", default="no")

if len(sys.argv)==1:
    parser.print_help()
    parser.exit()
(options, args) = parser.parse_args()

# CHECK input
if not os.path.isdir(options.path):
    sys.stderr.write("Found no directory with that name, {}\n".format(options.path))
    sys.exit()

# Create list of maf files.
maf_l = []
for fn in os.listdir(options.path):
	if fn.startswith(options.prefix):
		maf_l.append(fn)
maf_l = sorted(maf_l)

if bool(maf_l) == False:
	sys.stderr.write("Found no files with the given prefix, {}\n".format(options.prefix))
	sys.exit()

# Checking if the path ends with '/'
if (not options.path.endswith('/')):
	options.path = options.path+'/'

## if path not exist makedir pruned else not do it
if not os.path.isdir(options.removed):
    os.system('mkdir ' + str(options.removed))

# Iterate over  list.
for file_name in maf_l:

	print("Processing file: {}".format(file_name))

	# Open maf file.
	maf_input = open(options.path + file_name, 'r')

	# Create output files
	outfile1 = open(options.removed + '/rmOS_' + file_name, 'w')

	# Iterate over lines in file.
	add_order = 0
	for line in maf_input:
        line = line.strip()
		if line.startswith('#'):
			outfile1.write(line)
			add_order = 0
		elif line.startswith('s'):
			if '+' in line:
				if add_order == 0:
					outfile1.write('a ordered=true\n')
					add_order = 1
				outfile1.write(line)
		elif line == '\n':
			outfile1.write(line)
			add_order = 0

if options.clean=='yes':
	os.system('rm -r ' + str(options.path))

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_remove_opposite_strand.txt', 'x')
indication.close()
os.rename('./finished_remove_opposite_strand.txt', './output/finished_remove_opposite_strand.txt')
