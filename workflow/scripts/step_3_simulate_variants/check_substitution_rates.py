#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
:Author: Seyan Hu
:Date: 4-11-2022
:Extension and modification: Julia Höglund
:Date: 17-5-2023

:Usage: python <script.py> -i <name and path sim. variants> -l <path to logfiles>

To check for the nt substitution rate in the simulated variants file (vcf).
Validation step to determine whenever the simulation is performed correctly.  

'''

# Import dependencies.
import sys, os
from argparse import ArgumentParser
import natsort

parser = ArgumentParser(description=__doc__)
parser.add_argument("--sim-snps",
    help = "input file (genome-wide) with simulated variants ",
    type = str,
    required = True)
parser.add_argument("--trimmed-snps",
    help = "input file with (genome-wide) simulated variants filtered for ancestral overlap",
    type = str,
    required = True)
parser.add_argument("--param-logfiles",
    help = "input parameter logfiles created with 'process_parameters.py'", 
    type = str, 
    required = True,
    nargs="+")
parser.add_argument("--snp-outfile",
    help = "output file with snp summary", 
    type = str, 
    required = True)
parser.add_argument("--trimmed-outfile",
    help = "output file with filtered snp summary", 
    type = str, 
    required = True)
parser.add_argument("--param-outfile",
    help = "output parameter file with parameter summary", 
    type = str, 
    required = True)


args = parser.parse_args()

# Count for mutations. 
mut = 0
mutCpG = 0
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

## fill for later
chrom_list = []

#####################################
####### generated parameters ########
#####################################

def process_logfile(file, chrom_list):
    with open(file) as log_open:
        chrom = file.split('.')[0].split('chr')[1]
        chrom_list.append(chrom)

        for line in log_open:
            if line.startswith("#y\t"):
                line = log_open.readline()
                new_line = line.strip().split("\t")

                mut = int(new_line[0])
                ACn = int(new_line[2]) 
                AGn = int(new_line[3]) 
                ATn = int(new_line[4]) 

                CAn = int(new_line[5]) 
                CGn = int(new_line[6]) 
                CTn = int(new_line[7]) 

                GAn = int(new_line[8]) 
                GCn = int(new_line[9]) 
                GTn = int(new_line[10]) 

                TAn = int(new_line[11]) 
                TCn = int(new_line[12]) 
                TGn = int(new_line[13]) 

            elif line.startswith("#yCpG"):
                line = log_open.readline()
                new_line = line.strip().split("\t")
                mutCpG = int(new_line[0])

                CA = int(new_line[2]) 
                CG = int(new_line[3]) 
                CT = int(new_line[4]) 

                GA = int(new_line[5]) 
                GC = int(new_line[6]) 
                GT = int(new_line[7]) 

    return chrom, mut, ACn, AGn, ATn, CAn, CGn, CTn, GAn, GCn, GTn, TAn, TCn, TGn, mutCpG, CA, CG, CT, GA, GC, GT

def write_summary(log_file, chrom, mut, ACn, AGn, ATn, CAn, CGn, CTn, GAn, GCn, GTn, TAn, TCn, TGn, mutCpG, CA, CG, CT, GA, GC, GT):
    log_file.write("\nChromosome considered: " + str(chrom) + '\n')
    log_file.write('Total mut: ' + str(mut) + '\n')
    log_file.write("#AC\tAG\tAT\n")
    log_file.write(str(ACn) + '\t' + str(AGn) + '\t' + str(ATn) + '\n')
    log_file.write(str((100/mut)*ACn)+'%\t' + str((100/mut)*AGn)+'%\t' + str((100/mut)*ATn) + '%\n')
    log_file.write('#CA\tCG\tCT\n')
    log_file.write(str(CAn) + '\t' + str(CGn) + '\t' + str(CTn) + '\n')
    log_file.write(str((100/mut)*CAn) + '%\t' + str((100/mut)*CGn) + '%\t' + str((100/mut)*CTn) + '%\n')
    log_file.write('#GA\tGC\tGT\n')
    log_file.write(str(GAn) + '\t' + str(GCn) + '\t' + str(GTn) + '\n')
    log_file.write(str((100/mut)*GAn) + '%\t' + str((100/mut)*GCn) + '%\t' + str((100/mut)*GTn) + '%\n')
    log_file.write('#TA\tTC\tTG\n')
    log_file.write(str(TAn) + '\t' + str(TCn) + '\t' + str(TGn) + '\n')
    log_file.write(str((100/mut)*TAn) + '%\t' + str((100/mut)*TCn) + '%\t' + str((100/mut)*TGn) + '%\n')
    log_file.write('Total CpG mut: ' + str(mutCpG) + '\n')
    log_file.write("#CA\tCG\tCT\t(CpG)\n")
    log_file.write(str(CA) + '\t' + str(CG) + '\t' + str(CT) + '\n')
    log_file.write(str((100/mutCpG)*CA) + '%\t' + str((100/mutCpG)*CG) + '%\t' + str((100/mutCpG)*CT) + '%\n')
    log_file.write("#GA\tGC\tGT\t(CpG)\n")
    log_file.write(str(GA) + '\t' + str(GC) + '\t' + str(GT) + '\n')
    log_file.write(str((100/mutCpG)*GA) + '%\t' + str((100/mutCpG)*GC) + '%\t' + str((100/mutCpG)*GT) + '%\n')

log_logfiles = open(args.param_outfile, "w")
filenames = natsort.natsorted(args.param_logfiles)

print('compiling summary of created simulation parameters')
for file in filenames:
    chrom, mut, ACn, AGn, ATn, CAn, CGn, CTn, GAn, GCn, GTn, TAn, TCn, TGn, mutCpG, CA, CG, CT, GA, GC, GT = process_logfile(file, chrom_list)
    write_summary(log_logfiles, chrom, mut, ACn, AGn, ATn, CAn, CGn, CTn, GAn, GCn, GTn, TAn, TCn, TGn, mutCpG, CA, CG, CT, GA, GC, GT)

log_logfiles.close()

#####################################
######## simulated variants #########
#####################################

## reset counters
# Count for mutations. 
mut,mutCpG = 0,0
# Total nt substitution counts.
ACn, AGn, ATn, CAn, CGn, CTn, GAn, GCn, GTn, TAn, TCn, TGn = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
# Total nt substitution counts in CpG sites.
CA, CG, CT, GA, GC, GT = 0, 0, 0, 0, 0, 0

def process_vcf(vcf_file, chrom_list, log_file):
    for c in chrom_list:
        write = 0
        for line in vcf_file:
            if not line.startswith('#'):
                edited_line = line.strip().split('\t')
                chrom = edited_line[0]
                ref = edited_line[3]
                alt = edited_line[4]
                CpG_site = edited_line[7]

                if chrom == c:
                    # If there is no CpG site.
                    if not CpG_site == 'CpG':

                        if ref == 'A':
                            if alt == 'C':
                                ACn += 1
                                mut += 1
                            elif alt == 'G':
                                AGn += 1
                                mut += 1
                            elif alt == 'T':
                                ATn += 1
                                mut += 1

                        elif ref == 'C':
                            if alt == 'A':
                                CAn += 1
                                mut += 1
                            elif alt == 'G':
                                CGn += 1
                                mut += 1
                            elif alt == 'T':
                                CTn += 1
                                mut += 1
                    
                        elif ref == 'G':
                            if alt == 'A':
                                GAn += 1
                                mut += 1
                            elif alt == 'C':
                                GCn += 1
                                mut += 1
                            elif alt == 'T':
                                GTn += 1
                                mut += 1
                        
                        elif ref == 'T':
                            if alt == 'A':
                                TAn += 1
                                mut += 1
                            elif alt == 'C':
                                TCn += 1
                                mut += 1
                            elif alt == 'G':
                                TGn += 1
                                mut += 1
                    
                    # If there is a CpG site.
                    elif CpG_site == 'CpG':
                    
                        if ref == 'C':
                            if alt == 'A':
                                CA +=1
                                mutCpG += 1
                            elif alt == 'G':
                                CG += 1
                                mutCpG += 1
                            elif alt == 'T':
                                CT += 1
                                mutCpG += 1
                    
                        elif ref == 'G':
                            if alt == 'A':
                                GA += 1
                                mutCpG += 1
                            elif alt == 'C':
                                GC += 1
                                mutCpG += 1
                            elif alt == 'T':
                                GT += 1
                                mutCpG += 1
                    write = 1

        if vcf_file.readline() == '' and write == 1:
            ## Print data. 
            log_file.write("\nChromosome considered: " + str(c) + '\n')
            # For non CpG site mutations. 
            log_file.write('Total mut: ' + str(mut) + '\n')
            log_file.write("#AC\tAG\tAT\n")
            log_file.write(str(ACn) + '\t' + str(AGn) + '\t' + str(ATn) + '\n')
            log_file.write(str((100/mut)*ACn) + '%\t' + str((100/mut)*AGn) + '%\t' + str((100/mut)*ATn) + '%\n')
            
            log_file.write('#CA\tCG\tCT\n')
            log_file.write(str(CAn) + '\t' + str(CGn) + '\t' + str(CTn) + '\n')
            log_file.write(str((100/mut)*CAn) + '%\t'+ str((100/mut)*CGn) + '%\t' + str((100/mut)*CTn) + '%\n')
                
            log_file.write('#GA\tGC\tGT\n')
            log_file.write(str(GAn) + '\t' + str(GCn) + '\t' + str(GTn) + '\n')
            log_file.write(str((100/mut)*GAn) + '%\t' + str((100/mut)*GCn) + '%\t' + str((100/mut)*GTn) + '%\n')
                
            log_file.write('#TA\tTC\tTG\n')
            log_file.write(str(TAn) + '\t' + str(TCn) + '\t' + str(TGn) + '\n')
            log_file.write(str((100/mut)*TAn) + '%\t' + str((100/mut)*TCn) + '%\t' + str((100/mut)*TGn) + '%\n')
                
            # For CpG site mutations
            log_file.write('Total CpG mut: ' + str(mutCpG) + '\n')
                
            log_file.write("#CA\tCG\tCT\t(CpG)\n")
            log_file.write(str(CA) + '\t' + str(CG) + '\t' + str(CT) + '\n')
            log_file.write(str((100/mutCpG)*CA)+ '%\t' + str((100/mutCpG)*CG) + '%\t' + str((100/mutCpG)*CT) + '%\n')
                
            log_file.write("#GA\tGC\tGT\t(CpG)\n")
            log_file.write(str(GA) + '\t' + str(GC) + '\t' + str(GT) + '\n')
            log_file.write(str((100/mutCpG)*GA) + '%\t' + str((100/mutCpG)*GC) + '%\t' + str((100/mutCpG)*GT) + '%\n')  
        
        vcf_file.seek(0)

with open(args.sim_snps, 'r') as vcf_open, open(args.snp_outfile, "w") as log_simulated:
    print('compiling summary of simulated variants')
    process_vcf(vcf_open, chrom_list, log_simulated)

##########################################
######## filtered simulated snps #########
##########################################

## reset counters
# Count for mutations. 
mut,mutCpG = 0,0
# Total nt substitution counts.
ACn, AGn, ATn, CAn, CGn, CTn, GAn, GCn, GTn, TAn, TCn, TGn = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
# Total nt substitution counts in CpG sites.
CA, CG, CT, GA, GC, GT = 0, 0, 0, 0, 0, 0

with open(args.trimmed_snps, 'r') as filtered_open, open(args.trimmed_outfile, "w") as log_filtered:
    print('compiling summary of simulated variants filter for ancestor sequence coverage')
    process_vcf(filtered_open, chrom_list, log_filtered)
