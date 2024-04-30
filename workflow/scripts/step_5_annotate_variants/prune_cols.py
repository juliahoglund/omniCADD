#!/usr/bin/env python
"""
Prune certain columns from alignment
"""


import sys
import Bio
from Bio import AlignIO
                                   
__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"

def prune_aln(ofile, fh_out):
    """Prune what columns from alignment and print result
    """
    keep_cols = []
    try:
        for aln in AlignIO.parse(ofile, "fasta"):
            for i in range(aln.get_alignment_length()):
                # deprecated: col = aln.get_column(i)
                col_nucs = [sr.seq[i].upper() for sr in aln]
                if aln[0].seq[i]!='-':
                    keep_cols.append(i)
            for s in aln:
                fh_out.write(">%s\n" %(s.id))
                fh_out.write('%s\n' % ''.join([s.seq[i] for i in keep_cols]))
    except:
        pass


def main():
    """
    The main function
    """
    mfile = sys.argv[1]  # fasta file
    ofile = gzip.open(mfile, "rt") \
                if mfile.endswith('.gz') else open(mfile, "r") # open file
    file_h = sys.argv[2]
    outfile = open(file_h, 'w')
    prune_aln(ofile, outfile)


if __name__ == "__main__":
    main()

# alignments = list(AlignIO.parse(ofile, "maf"))