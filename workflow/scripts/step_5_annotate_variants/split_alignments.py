#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Split a multi-FASTA alignment file into smaller chunks.
Derived from: split_fasta.py (from: https://thevirtuallaboratory.com/blog/splitting-a-multi-fasta)
Author: Bram van Dijk
Modified by: Julia Höglund 2025-03-05
"""

import sys
import math
import gzip
import Bio
from Bio import AlignIO
import os

def batch_iterator(iterator, batch_size):
    """
    This is a generator function, and it returns lists of the
    entries from the supplied iterator. Each list will have
    batch_size entries, although the final list may be shorter.
    (derived from https://biopython.org/wiki/Split_large_file)
    """
    entry = True  # Make sure we loop once
    iter_object = iter(iterator)
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iter_object)
            except StopIteration:
                entry = None
            if entry is None:
                break  # EOF = end of file
            batch.append(entry)
        if batch:
            yield batch

if len(sys.argv) != 5:
    sys.exit("usage: convert_alignments.py MAF_FILE N_CHUNKS OUTPUT_FOLDER_MAF REF_SPECIES")

mfile = sys.argv[1]  # maf file
ofile = gzip.open(mfile, "rt") if mfile.endswith('.gz') else open(mfile, "r")  # maf file

chunks = sys.argv[2]  # number of chunks
maf_folder = sys.argv[3]  # folder to save maf chunks in
ref_species = sys.argv[4]

if not os.path.exists(maf_folder):
    os.makedirs(maf_folder)

to_keep = []

# parse alignment
for alignment in AlignIO.parse(ofile, "maf"):
    if ref_species in str(alignment):
        to_keep.append(alignment)

nseq = len(to_keep)
chunksize = math.ceil(nseq / int(chunks))
chrom = mfile.split('.')[0].split('chr')[1]

print("Splitting maffile file of", nseq, "blocks into chunks of", chunksize, "blocks")
for i, batch in enumerate(batch_iterator(to_keep, chunksize)):
    filename = os.path.join(maf_folder, f"chr{chrom}-{i + 1}.maf")
    with open(filename, "w") as maf_handle:
        count = Bio.AlignIO.write(batch, maf_handle, "maf")
    print("Wrote %i sequences to %s" % (count, filename))
