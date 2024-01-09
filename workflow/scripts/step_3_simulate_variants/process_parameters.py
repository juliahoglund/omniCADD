#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

# :Author: Martin Kircher
# :Contact: mkircher@uw.edu
# :Date: *09.05.2012
# :Reformatted by: Julia Höglund
# :Contact: julia.hoglund@su.se
# :Date: *21.11.2022
# :Reformatted by: Job van Schipstal
# :Contact: job.vanschipstal@wur.nl
# :Date: *27.09.2023

:Usage: see process_parameters.py --help

This script processes the logfiles produced by create_parameters.py so that
they are ready to use for simulate.variants.py

This functionality was split off from the simulate_variants script.
"""

import sys
from collections import defaultdict  # make dicts that can be filled on the go
from argparse import ArgumentParser
import traceback
import pickle

parser = ArgumentParser(description=__doc__)
parser.add_argument("-n", "--number",
    help="Approximate number of substitution events to simulate (default: 1000)", 
    type=str, 
    default="1000")
parser.add_argument("-p", "--params",
    help="Input parameter log files, 1 or more files",
    type=str, 
    required=True, 
    nargs="+")
parser.add_argument("-l", "--logfile",
    help="logfile (default process_parameters.txt)", 
    type=str,
    default="process_parameters.txt")
parser.add_argument("-o", "--outfile",
    help="Output pickle file containing the needed simulation parameters for simulate_variants.py",
    type=str, 
    required=True)

args = parser.parse_args()

# CREATE COUNTERS OF STATS AND COMPUTATION
totalrefA = 0
totalrefC = 0
totalrefG = 0
totalrefT = 0

mut, total = 0, 0
nAC, nAG, nAT, nCA, nCG, nCT, nGA = 0, 0, 0, 0, 0, 0, 0
nGC, nGT, nTA, nTC, nTG = 0, 0, 0, 0, 0
totalCpG, mutCpG = 0, 0
mCA, mCG, mCT, mGA, mGC, mGT = 0, 0, 0, 0, 0, 0

insertionsizes, deletionsizes = defaultdict(int), defaultdict(int)

intervals = defaultdict(list)
interval_vals = {}

# FILL COUNTERS; COUNTS AND PRINTS STATS
control_mut, control_mutCpG, control_total = 0, 0, 0
try:
    for filename in args.params:
        infile = open(filename)
        sys.stderr.write("Logfile considered: %s\n" % filename)
        line = infile.readline()
        while line != "":
            if line == "#A\tC\tG\tT\tCpGs\n":
                line = infile.readline()
                A, C, G, T, cCpG = list(map(int, line.split()))
                totalrefA += A
                totalrefC += C
                totalrefG += G
                totalrefT += T
                totalCpG += cCpG
            elif line == "#y\tN\tAC\tAG\tAT\tCA\tCG\tCT\tGA\tGC\tGT\tTA\tTC\tTG\n":
                line = infile.readline()
                cmut, ctotal, cAC, cAG, cAT, cCA, cCG, cCT, cGA, cGC, \
                cGT, cTA, cTC, cTG = list(map(int, line.split()))
                mut += cmut
                total += ctotal
                nAC += cAC
                nAG += cAG
                nAT += cAT
                nCA += cCA
                nCG += cCG
                nCT += cCT
                nGA += cGA
                nGC += cGC
                nGT += cGT
                nTA += cTA
                nTC += cTC
                nTG += cTG
            elif line == "#yCpG\tNCpG\tCA\tCG\tCT\tGA\tGC\tGT\n":
                line = infile.readline()
                cmutCpG, ctotalCpG, cCA, cCG, cCT, cGA, cGC, cGT = list(map(int, line.split()))
                mutCpG += cmutCpG
                totalCpG += ctotalCpG
                mCA += cCA
                mCG += cCG
                mCT += cCT
                mGA += cGA
                mGC += cGC
                mGT += cGT
            elif line == "##INSERTIONS\n":
                infile.readline()  # HEADER
                line = infile.readline()
                while line != "" and not line.startswith('#'):
                    fields = line.split()
                    if len(fields) == 2:
                        length, count = list(map(int, fields[:2]))
                        insertionsizes[length] += count
                    line = infile.readline()
            elif line == "##DELETIONS\n":
                infile.readline()  # HEADER
                line = infile.readline()
                while line != "" and not line.startswith('#'):
                    fields = line.split()
                    if len(fields) == 2:
                        length, count = list(map(int, fields[:2]))
                        deletionsizes[length] += count
                    line = infile.readline()
            else:
                fields = line.split()
                if len(fields) == 11:
                    chrom = fields[0]
                    window_start, window_end, cmut, ctotal, cmutCpG, ctotalCpG, cA, cC, cG, cT = list(map(int, fields[1:11]))
                    control_mut += cmut
                    control_mutCpG += cmutCpG
                    control_total += ctotal
                    intervals[chrom].append((window_start, window_end + 1))
                    interval_vals[(chrom, window_start, window_end + 1)] = (cmut, ctotal, cmutCpG, ctotalCpG, cA, cC, cG, cT)
                line = infile.readline()
        infile.close()
except:
    exc_type, exc_value, exc_traceback = sys.exc_info()
    sys.stderr.write("%s\n" % str(exc_value))
    traceback.print_tb(exc_traceback)
    sys.stderr.write('Script terminated early. Printing current values.\n')

# PRINT STATS TO LOGFILE
logfile = open(args.logfile, 'w')

logfile.write("Chromosomes considered: %s\n" % len(intervals))
logfile.write("Number of regions: %d\n" % len(interval_vals))

sum_obs = float(totalrefA + totalrefC + totalrefG + totalrefT)
total = float(total)
totalCpG = float(totalCpG)
logfile.write("Total aligned positions: %d (Exp: %d)\n" % (sum_obs, control_total))
logfile.write("Total mutations: %d/%d (Exp: %d/%d); Rate: %.8f/%.8f\n" % (
    mut, mutCpG, control_mut, control_mutCpG, mut / total, mutCpG / totalCpG))
gmut = mut / total
gmutCpG = mutCpG / totalCpG

gmut_ = mut / sum_obs
gmutCpG_ = mutCpG / sum_obs
logfile.write(
    "Abs. mutation rate non-CpG: %.8f\tAbs. mutation rate CpGs: %.8f\n" % (
        gmut_, gmutCpG_))

total_variants = int(float(args.number))  # To support scientific notation
dmut = (total_variants * gmut_ / (gmut_ + gmutCpG_)) / total
dmutCpG = (total_variants * gmutCpG_ / (gmut_ + gmutCpG_)) / totalCpG
dmutindel = total_variants / sum_obs
logfile.write("Likelihood altering a CpG: %.8f; Altering another base: %.8f\n" % (dmutCpG, dmut))

# RATES DETERMINED FROM HUMAN-CHIMP ANCESTOR (CORRECTED FOR BASE COMPOSITION)

totalrefA = float(totalrefA)
totalrefC = float(totalrefC)
totalrefG = float(totalrefG)
totalrefT = float(totalrefT)

# SUBSTITUTION MATRIX NON-CpG
GTR_nonCpG = {'A': [              0, nAC / totalrefA, nAG / totalrefA, nAT / totalrefA],
              'C': [nCA / totalrefC,               0, nCG / totalrefC, nCT / totalrefC],
              'G': [nGA / totalrefG, nGC / totalrefG,               0, nGT / totalrefG],
              'T': [nTA / totalrefT, nTC / totalrefT, nTG / totalrefT,               0]}

ffactor = [1.0 / sum(x) for x in list(GTR_nonCpG.values())]
for ki, key in enumerate(GTR_nonCpG.keys()):
    for i in range(4):
        GTR_nonCpG[key][i] = GTR_nonCpG[key][i] * ffactor[ki]

# SUBSTITUTION MATRIX CpG
mtotalrefC = (totalCpG - (
        mCA + mCG + mCT + mGA + mGC + mGT)) / 2.0 + mCA + mCG + mCT
mtotalrefG = (totalCpG - (
        mCA + mCG + mCT + mGA + mGC + mGT)) / 2.0 + mGA + mGC + mGT

GTR_CpG = {'C': [mCA / mtotalrefC,                0, mCG / mtotalrefC, mCT / mtotalrefC],
           'G': [mGA / mtotalrefG, mGC / mtotalrefG,                0, mGT / mtotalrefG]}

ffactor = [1.0 / sum(x) for x in list(GTR_CpG.values())]
for ki, key in enumerate(GTR_CpG.keys()):
    for i in range(4):
        GTR_CpG[key][i] = GTR_CpG[key][i] * ffactor[ki]

# INITIATE INSERT NORMALIZED LIKELIHOODS
inserts = []
tinserts = float(sum(insertionsizes.values()))
logfile.write("TOTAL INSERTS FROM MSA:\t%d\n" % tinserts)
to_sort = list(insertionsizes.keys())
to_sort.sort()
logfile.write("Length\tCount\tnFreq\n")
for key in to_sort:
    value = insertionsizes[key] / tinserts
    inserts.append((key, value))
del insertionsizes
ginserts = tinserts / sum_obs
logfile.write(
    "INSERTIONS: gfreq=%.8f, %s...\n" % (ginserts, str(inserts)[:120]))

# INITIATE DELETION NORMALIZED LIKELIHOODS
deletions = []
tdeletions = float(sum(deletionsizes.values()))
logfile.write("TOTAL DELETIONS FROM MSA:\t%d\n" % tdeletions)
to_sort = list(deletionsizes.keys())
to_sort.sort()
logfile.write("Length\tCount\tnFreq\n")
for key in to_sort:
    value = deletionsizes[key] / tdeletions
    # if key <= 10: logfile.write("%d\t%d\t%.8f\n"%(key, deletionsizes[key],
    # value))
    deletions.append((key, value))
del deletionsizes
gdeletions = tdeletions / sum_obs
logfile.write(
    "DELETIONS: gfreq=%.8f, %s...\n" % (gdeletions, str(deletions)[:120]))

pickle_f = open(args.outfile, 'wb')
pickle.dump([intervals, interval_vals, gmut, dmut, gmutCpG, dmutCpG,
             gdeletions, ginserts, dmutindel, GTR_CpG, GTR_nonCpG,
             deletions, inserts], pickle_f)
pickle_f.close()
