#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

# :Author: Martin Kircher
# :Contact: mkircher@uw.edu
# :Date: *09.05.2012
# :Reformatted by: Julia Höglund
# :Contact: julia.hoglund@su.se
# :Date: *21.11.2022

:Usage: python <script.py> -p <folder containing log files> -i <infile prefix of genome> -o <name of outfile> -c <chromosomes> -n <number of events>

:Example:
python scripts/apply_parameters.py -p ./output/logfiles -i ./genome/sus_scrofa_ref.fa -o sim_variants.vcf -c 1,2,3,4,5 -n 10000

This script simulates the variants, based on the log files that has been created with create_parameters.py.

"""

import os,sys
from collections import defaultdict # make dicts that can be filled on the go
from optparse import OptionParser   # make option parser
import traceback
import random                       # create random numbers from distr.
import mmap                         # allow for memory mapping
from progressbar import Percentage, ProgressBar, Bar, ETA # allow progress bar

# READ FASTA INDEX TO MEMORY
def read_fasta_index(filename):
  infile = open(filename)
  res = {}
  for line in infile:
    fields = line.split()
    if len(fields) == 5:
      cname,length,start,line,cline = fields
      res[cname]=int(length),int(start),int(line),int(cline)
    else:
      sys.stderr.write('Error: Unexpected line in fasta index file: %s\n'%(line.strip()))
      sys.exit()
  infile.close()
  return res

# CREATE OPTION PARSER
parser = OptionParser()
parser.add_option("-n","--number_events", dest = "number", help = "Approximate number of substitution events to simulate (default 1000)",type = "int", default = 1000)
parser.add_option("-c","--chroms", dest = "chromosomes", help = "Chromosomes considered", default = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22')
parser.add_option("-p","--params_folder", dest = "params_folder", help = "Folder with parameter files and prefix of files (def output/logfiles/chr)", default = "output/logfiles")
parser.add_option("-i", "--infile", dest = "infile", help = "Infile name of genome (prefix for fasta index, def reference_genome.fa)", default = "reference_genome.fa")
parser.add_option("-o", "--outfile", dest = "outfile", help = "Name of output file (def STDOUT)", default = None)
(options, args) = parser.parse_args()

# OPEN OUTPUT FILE IF GIVEN, OTHERWISE PRINT ON SCREEN
if options.outfile != None:
  output = open(options.outfile,'w')
else:
  output = sys.stdout

# PARSE GENOME FILE(s)
if not os.path.exists(options.infile) or not os.path.exists(options.infile+".fai"):
  sys.stderr.write("Invalid path for genome and genome fasta index file.\n")
  sys.exit()

fastaindex = read_fasta_index(options.infile+".fai")

# PARSE MSA LOGFILE(s)
if not os.path.isdir('/'.join(options.params_folder.split('/')[:-1])):
  sys.stderr.write("Invalid path for parameter (log) files.\n")
  sys.exit()

if options.chromosomes != '':
  filenames = [options.params_folder+"%s.log"%x for x in options.chromosomes.split(',')]
else:
  filenames = [x for x in os.path.listdir(options.params_folder) if x.endswith(".log")]

# CREATE COUNTERS OF STATS AND COMPUTATION
totalrefA = 0
totalrefC = 0
totalrefG = 0
totalrefT = 0

mut,total = 0,0
nAC,nAG,nAT,nCA,nCG,nCT,nGA,nGC,nGT,nTA,nTC,nTG = 0,0,0,0,0,0,0,0,0,0,0,0
totalCpG,mutCpG = 0,0
mCA,mCG,mCT,mGA,mGC,mGT = 0,0,0,0,0,0

insertionsizes,deletionsizes = defaultdict(int),defaultdict(int)

intervals = defaultdict(list)
interval_vals = {}

# FILL COUNTERS; COUNTS AND PRINTS STATS
control_mut,control_mutCpG,control_total = 0,0,0
try:
  for filename in filenames:
    infile = open(filename)
    sys.stderr.write("Logfile considered: %s\n"%filename)
    line = infile.readline()
    while line != "":
      if line == "#A\tC\tG\tT\tCpGs\n":
        line = infile.readline()
        A,C,G,T,cCpG = list(map(int,line.split()))
        totalrefA += A
        totalrefC += C
        totalrefG += G
        totalrefT += T
        totalCpG += cCpG
      elif line == "#y\tN\tAC\tAG\tAT\tCA\tCG\tCT\tGA\tGC\tGT\tTA\tTC\tTG\n":
        line = infile.readline()
        cmut,ctotal,cAC,cAG,cAT,cCA,cCG,cCT,cGA,cGC,cGT,cTA,cTC,cTG = list(map(int,line.split()))
        mut+=cmut
        total+=ctotal
        nAC+=cAC
        nAG+=cAG
        nAT+=cAT
        nCA+=cCA
        nCG+=cCG
        nCT+=cCT
        nGA+=cGA
        nGC+=cGC
        nGT+=cGT
        nTA+=cTA
        nTC+=cTC
        nTG+=cTG
      elif line == "#yCpG\tNCpG\tCA\tCG\tCT\tGA\tGC\tGT\n":
        line = infile.readline()
        cmutCpG,ctotalCpG,cCA,cCG,cCT,cGA,cGC,cGT = list(map(int,line.split()))
        mutCpG+=cmutCpG
        totalCpG+=ctotalCpG
        mCA+=cCA
        mCG+=cCG
        mCT+=cCT
        mGA+=cGA
        mGC+=cGC
        mGT+=cGT
      elif line == "##INSERTIONS\n":
        line = infile.readline() # HEADER
        line = infile.readline()
        while line != "" and not line.startswith('#'):
          fields = line.split()
          if len(fields) == 2:
            length,count=list(map(int,fields[:2]))
            insertionsizes[length]+=count
          line = infile.readline()
      elif line == "##DELETIONS\n":
        line = infile.readline() # HEADER
        line = infile.readline()
        while line != "" and not line.startswith('#'):
          fields = line.split()
          if len(fields) == 2:
            length,count=list(map(int,fields[:2]))
            deletionsizes[length]+=count
          line = infile.readline()
      else:
        fields = line.split()
        if len(fields) == 11:
          chrom = fields[0]
          window_start,window_end,cmut,ctotal,cmutCpG,ctotalCpG,cA,cC,cG,cT = list(map(int,fields[1:11]))
          control_mut+=cmut
          control_mutCpG+=cmutCpG
          control_total+=ctotal
          intervals[chrom].append((window_start, window_end+1))
          interval_vals[(chrom,window_start, window_end+1)]=(cmut,ctotal,cmutCpG,ctotalCpG,cA,cC,cG,cT)
        line = infile.readline()
    infile.close()
except:
  exc_type, exc_value, exc_traceback = sys.exc_info()
  sys.stderr.write("%s\n"%str(exc_value))
  traceback.print_tb(exc_traceback)
  sys.stderr.write('Script terminated early. Printing current values.\n')

# PRINT STATS TO LOGFILE
logfile = open("parameter_logfile.log",'w')

logfile.write("Chromosomes considered: %s\n"%len(intervals))
logfile.write("Number of regions: %d\n"%len(interval_vals))

sum_obs = float(totalrefA+totalrefC+totalrefG+totalrefT)
total = float(total)
totalCpG = float(totalCpG)
logfile.write("Total aligned positions: %d (Exp: %d)\n"%(sum_obs,control_total))
logfile.write("Total mutations: %d/%d (Exp: %d/%d); Rate: %.8f/%.8f\n"%(mut,mutCpG,control_mut,control_mutCpG,mut/total,mutCpG/totalCpG))
gmut = mut/total
gmutCpG = mutCpG/totalCpG

gmut_ = mut/sum_obs
gmutCpG_ = mutCpG/sum_obs
logfile.write("Abs. mutation rate non-CpG: %.8f\tAbs. mutation rate CpGs: %.8f\n"%(gmut_,gmutCpG_))

dmut = (options.number*gmut_/(gmut_+gmutCpG_))/total
dmutCpG = (options.number*gmutCpG_/(gmut_+gmutCpG_))/totalCpG
dmutindel = options.number/sum_obs
logfile.write("Likelihood altering a CpG: %.8f; Altering another base: %.8f\n"%(dmutCpG,dmut))

# RATES DETERMINED FROM HUMAN-CHIMP ANCESTOR (CORRECTED FOR BASE COMPOSITION)

totalrefA = float(totalrefA)
totalrefC = float(totalrefC)
totalrefG = float(totalrefG)
totalrefT = float(totalrefT)

##SUBSTITUTION MATRIX NON-CpG
GTR_nonCpG = {'A':[             0 , nAC/totalrefA , nAG/totalrefA , nAT/totalrefA ],
              'C':[ nCA/totalrefC ,             0 , nCG/totalrefC , nCT/totalrefC ],
              'G':[ nGA/totalrefG , nGC/totalrefG ,             0 , nGT/totalrefG ],
              'T':[ nTA/totalrefT , nTC/totalrefT , nTG/totalrefT ,             0 ] }

ffactor = [1.0/sum(x) for x in list(GTR_nonCpG.values())]
for ki,key in enumerate(GTR_nonCpG.keys()):
  for i in range(4):
    GTR_nonCpG[key][i] = GTR_nonCpG[key][i] * ffactor[ki]

##SUBSTITUTION MATRIX CpG
mtotalrefC = (totalCpG-(mCA+mCG+mCT+mGA+mGC+mGT))/2.0+mCA+mCG+mCT
mtotalrefG = (totalCpG-(mCA+mCG+mCT+mGA+mGC+mGT))/2.0+mGA+mGC+mGT

GTR_CpG = {   'C':[ mCA/mtotalrefC ,              0 , mCG/mtotalrefC , mCT/mtotalrefC ],
              'G':[ mGA/mtotalrefG , mGC/mtotalrefG ,              0 , mGT/mtotalrefG ] }

ffactor = [1.0/sum(x) for x in list(GTR_CpG.values())]
for ki,key in enumerate(GTR_CpG.keys()):
  for i in range(4):
    GTR_CpG[key][i] = GTR_CpG[key][i] * ffactor[ki]

# INITIATE INSERT NORMALIZED LIKELIHOODS
inserts = []
tinserts = float(sum(insertionsizes.values()))
logfile.write("TOTAL INSERTS FROM MSA:\t%d\n"%tinserts)
to_sort = list(insertionsizes.keys())
to_sort.sort()
logfile.write("Length\tCount\tnFreq\n")
for key in to_sort:
  value = insertionsizes[key]/tinserts
  inserts.append((key,value))
del insertionsizes
ginserts = tinserts/sum_obs
logfile.write("INSERTIONS: gfreq=%.8f, %s...\n"%(ginserts,str(inserts)[:120]))

# INITIATE DELETION NORMALIZED LIKELIHOODS
deletions = []
tdeletions = float(sum(deletionsizes.values()))
logfile.write("TOTAL DELETIONS FROM MSA:\t%d\n"%tdeletions)
to_sort = list(deletionsizes.keys())
to_sort.sort()
logfile.write("Length\tCount\tnFreq\n")
for key in to_sort:
  value = deletionsizes[key]/tdeletions
  #if key <= 10: logfile.write("%d\t%d\t%.8f\n"%(key, deletionsizes[key], value))
  deletions.append((key,value))
del deletionsizes
gdeletions = tdeletions/sum_obs
logfile.write("DELETIONS: gfreq=%.8f, %s...\n"%(gdeletions,str(deletions)[:120]))

#OPEN FASTA COPY WITH MMAP
sys.stderr.write('Doing mmap of genome file...\n')
f = open(options.infile, "r")
cmap = mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_READ)
cmap.flush()

sys.stderr.write('Iterating over genome...\n')

# ITERATE; SIMULATE VARIANTS, PRINT TO VCF OUTPUT
VCFheader = """##fileformat=VCFv4.1
##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context (based on the reference sequence).
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"""
output.write(VCFheader+"\n")
to_sort = list(intervals.keys())
to_sort.sort()

for chrom in to_sort:
  if chrom in fastaindex:
    length,sblock,bline,cline = fastaindex[chrom]
    cintervals = intervals[chrom]
    cintervals.sort()
    sys.stderr.write("Simulating variants on chrom\t%s\n"%chrom)
    N = len(cintervals)
    pbar = ProgressBar(widgets=[Bar('=', '[', ']'), ' ', Percentage(), ' ', ETA()], maxval = N).start()
    def job():
       for ind,(start,end) in enumerate(cintervals):
         wmut,wmutCpG,wtotal,wtotalCpG,wA,wC,wG,wT = 0,0,0,0,0,0,0,0
         for neighbors in range( max(0,ind-5-max(ind+6-len(cintervals),0)), min(len(cintervals),ind+6-min(0,ind-5)) ):
           pos1,pos2 = cintervals[neighbors]
           (cmut,ctotal,cmutCpG,ctotalCpG,cA,cC,cG,cT) = interval_vals[(chrom,pos1,pos2)]
           wmut += cmut
           wtotal += ctotal
           wmutCpG += cmutCpG
           wtotalCpG += ctotalCpG
           wA += cA
           wC += cC
           wG += cG
           wT += cT
         if wtotal == 0:
           sys.stderr.write('Encountered genomic blocks without any non-N bases. Skipping.\n')
           continue
         tbases = float(wA+wC+wG+wT)
         bfreqs = [wA/tbases,wC/tbases,wG/tbases,wT/tbases]
         if wtotal > 1000: y = wmut/float(wtotal)
         else: y=gmut
         if wtotalCpG > 1000: yCpG = wmutCpG/float(wtotalCpG)
         else: yCpG=gmutCpG
         #print bfreqs,y,yCpG

         #dmut, gmut, gdeletions (deletions), ginserts (inserts)
         lmut = y/gmut * dmut
         lmutGpG = yCpG/gmutCpG * dmutCpG
         ldeletion = y/gmut * gdeletions/gmut * dmutindel
         linsertion = y/gmut * ginserts/gmut * dmutindel

         if not end+1 <= length: end = length - 1
         if not start > 1: start = 2

         lbase,base,nbase = None,None,None
         for pos in range(start,end):
           if lbase != None:
             lbase,base,nbase = base,nbase,None
             npos = (pos // bline)*cline+(pos % bline)
             cmap.seek(sblock+npos)
             nbase = cmap.read(1).decode("utf-8").upper()
           else:
             npos = ((pos-2) // bline)*cline+((pos-2) % bline)
             cmap.seek(sblock+npos)
             lbase = cmap.read(1).decode("utf-8").upper()
             npos = ((pos-1) // bline)*cline+((pos-1) % bline)
             cmap.seek(sblock+npos)
             base = cmap.read(1).decode("utf-8").upper()
             npos = (pos // bline)*cline+(pos % bline)
             cmap.seek(sblock+npos)
             nbase = cmap.read(1).decode("utf-8").upper()

           ismut = random.random()
           if (base == "C" and nbase == "G") or (base == "G" and lbase == "C"):
             if ismut <= lmutGpG:
               whichmut = random.random()
               for i,alt in enumerate('ACGT'):
                 if alt != base:
                   whichmut-=GTR_CpG[base][i]
                   if whichmut < 0:
                     output.write("%s\t%d\t.\t%s\t%s\t.\t.\tCpG\t.\n"%(chrom,pos,base,alt))
                     break
           elif base in "ACGT":
             if ismut <= lmut:
               whichmut = random.random()
               for i,alt in enumerate('ACGT'):
                 if alt != base:
                   whichmut-=GTR_nonCpG[base][i]
                   if whichmut < 0:
                     output.write("%s\t%d\t.\t%s\t%s\t.\t.\t.\t.\n"%(chrom,pos,base,alt))
                     break
           else:
             #sys.stderr.write("Base %s on %s %d can not be mutated. Skipping.\n"%(base,chrom,pos))
             continue

           isindel = random.random()
           if isindel <= ldeletion:
             npos = ((pos-1) // bline)*cline+((pos-1) % bline)
             cmap.seek(sblock+npos)
             base = cmap.read(1).decode("utf-8").upper()
             if base == "N": continue

             indellength = random.random()
             for dellength,value in deletions:
               indellength-=value
               if indellength <= 0:
                 npos = (pos // bline)*cline+(pos % bline)
                 cmap.seek(sblock+npos)
                 delseq = ''
                 for i in range(dellength):
                   if pos+i < length:
                     npos = ((pos+i) // bline)*cline+((pos+i) % bline)
                     cmap.seek(sblock+npos)
                     delseq += cmap.read(1).decode("utf-8").upper()
                   else:
                     break
                 if len(delseq) == dellength:
                   output.write("%s\t%d\t.\t%s\t%s\t.\t.\t.\t.\n"%(chrom,pos,base+delseq,base))
                 break

           isindel = random.random()
           if isindel <= linsertion:
             npos = ((pos-1) // bline)*cline+((pos-1) % bline)
             cmap.seek(sblock+npos)
             base = cmap.read(1).decode("utf-8").upper()
             if base == "N": continue

             indellength = random.random()
             for inslength,value in inserts:
               indellength-=value
               if indellength <= 0:
                 insseq = ""
                 for s in range(inslength):
                   whichbase = random.random()
                   for i,newbase in enumerate('ACGT'):
                     whichbase -= bfreqs[i]
                     if whichbase <= 0:
                       insseq += newbase
                       break
                 output.write("%s\t%d\t.\t%s\t%s\t.\t.\t.\t.\n"%(chrom,pos,base,base+insseq))
                 break
         pbar.update(ind)
       pbar.finish()
       sys.stderr.write("Chrom %s done.\n"%chrom)
  else:
      sys.stderr.write('Chromosome %s not available from fasta index. Skipping.\n')
  job()
cmap.close()
f.close()

print("Step 2; Simulate variants done.\n")

# Create a txt file indicating that this process is finished (for snakemake)
indication = open('finished_apply_parameters.txt', 'x')
indication.close()
os.rename('./finished_apply_parameters.txt', './output/finished_apply_parameters.txt')
