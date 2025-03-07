#!/usr/bin/env python
"""
Prune certain columns from alignment
Modified by: Julia Höglund 2025-03-05
"""

import sys
import Bio
from Bio import AlignIO
import gzip

__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@gmail.com"
__license__ = "The MIT License (MIT)"

def prune_aln(ofile, fh_out):
    """Prune columns from alignment and print result."""
    keep_cols = []
    try:
        for aln in AlignIO.parse(ofile, "fasta"):
            for i in range(aln.get_alignment_length()):
                col_nucs = [sr.seq[i].upper() for sr in aln]
                if aln[0].seq[i] != '-':
                    keep_cols.append(i)
            for s in aln:
                fh_out.write(">%s\n" % (s.id))
                fh_out.write('%s\n' % ''.join([s.seq[i] for i in keep_cols]))
    except Exception as e:
        print(f"Error processing alignment: {e}")

def main():
    """Main function to handle file input and output."""
    mfile = sys.argv[1]  # fasta file
    file_h = sys.argv[2]  # output file

    with (gzip.open(mfile, "rt") if mfile.endswith('.gz') else open(mfile, "r")) as ofile, open(file_h, 'w') as outfile:
        prune_aln(ofile, outfile)

if __name__ == "__main__":
    main()