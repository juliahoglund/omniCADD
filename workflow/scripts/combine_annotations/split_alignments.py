#!/usr/bin/env python3
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
import lz4.frame
import Bio
from Bio import AlignIO
import os


def open_maf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith(".lz4"):
        return lz4.frame.open(path, "rt")
    return open(path, "r")


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
    sys.exit("usage: split_alignments.py MAF_FILE N_CHUNKS OUTPUT_FOLDER_MAF REF_SPECIES")

mfile = sys.argv[1]  # maf file
chunks = sys.argv[2]  # number of chunks
maf_folder = sys.argv[3]  # folder to save maf chunks in
ref_species = sys.argv[4]

if not os.path.exists(maf_folder):
    os.makedirs(maf_folder)

# First pass: count matching blocks without holding them in memory
nseq = 0
with open_maf(mfile) as handle:
    for alignment in AlignIO.parse(handle, "maf"):
        if ref_species in str(alignment):
            nseq += 1

chunksize = math.ceil(nseq / int(chunks))
chrom = mfile.split('.')[0].split('chr')[1]

print("Splitting maffile file of", nseq, "blocks into chunks of", chunksize, "blocks")

# Second pass: stream matching blocks straight into batches and write each
# batch out immediately, rather than accumulating the whole chromosome in
# memory before writing anything
with open_maf(mfile) as handle:
    matching = (alignment for alignment in AlignIO.parse(handle, "maf") if ref_species in str(alignment))
    for i, batch in enumerate(batch_iterator(matching, chunksize)):
        filename = os.path.join(maf_folder, f"chr{chrom}-{i + 1}.maf.lz4")
        with lz4.frame.open(filename, "wt") as maf_handle:
            count = Bio.AlignIO.write(batch, maf_handle, "maf")
        print("Wrote %i sequences to %s" % (count, filename))
