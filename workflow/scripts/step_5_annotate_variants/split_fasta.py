#!/usr/bin/env python
# -*- coding: utf-8 -*-

# split_fasta.py (assumes you have biopython installed, e.g. with pip install biopython)
# from: https://thevirtuallaboratory.com/blog/splitting-a-multi-fasta
# Author: Bram van Dijk

import sys, math
from Bio import SeqIO

def batch_iterator(iterator, batch_size):
    """
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    (derived from https://biopython.org/wiki/Split_large_file)
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                break # EOF = end of file
            batch.append(entry)
        if batch:
            yield batch

if(len(sys.argv) != 4):
        sys.exit("usage: split_fasta.py MULTI_FASTA_FILE N_CHUNKS OUTPUT_FOLDER")

ffile=sys.argv[1]  # fasta file
chunks=sys.argv[2] # number of chunks
folder=sys.argv[3] # folder to save chunks in

nseq = len([1 for line in open(ffile) if line.startswith(">")])
chunksize=math.ceil(nseq/int(chunks))
chrom = ffile.split('.')[1].split('chr')[2]
print("Splitting multi-fasta file of", nseq, "sequences into chunks of size", chunksize)

records = SeqIO.parse(open(ffile), "fasta")
for i, batch in enumerate(batch_iterator(records, chunksize)):
        filename = str(folder) + "chr" + chrom + "_%i.fasta" % (i + 1)
        with open(filename, "w") as handle:
                count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i sequences to %s" % (count, filename))

sys.exit("Done.")