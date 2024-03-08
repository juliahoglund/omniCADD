#!/usr/bin/env python
# -*- coding: utf-8 -*-

# echo sus_scrofa > species.list
# seqtk subseq results/alignment/fasta/43_mammals.epo/chr1_11.fasta species.list > test.fa
# grep -v ">" test.fa | fold -w1 > chr10_11.fa

# dependencies
import sys, math
import Bio
from Bio import AlignIO

def batch_iterator(iterator, batch_size):
    """
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
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
                break # EOF = end of file
            batch.append(entry)
        if batch:
            yield batch


if(len(sys.argv) != 5):
        sys.exit("usage: split_fasta.py MULTI_FASTA_FILE N_CHUNKS OUTPUT_FOLDER_MAF OUTPUT_FOLDER_FASTA")

ffile=sys.argv[1]  # maf file
chunks=sys.argv[2] # number of chunks
maf_folder=sys.argv[3] # folder to save chunks in
fasta_folder=sys.argv[4] # folder to save converted chunks in

to_keep = []

# parse alignment
for alignment in AlignIO.parse(ffile, "maf"):
    print("Alignment of length %i" % alignment.get_alignment_length())
    if "sus_scrofa" in str(alignment): 
        print(alignment)
        to_keep.append(alignment)

nseq = len(to_keep)
chunksize=math.ceil(nseq/int(chunks))
start = 0
chrom = ffile.split('.')[1].split('chr')[2]

print("Splitting maffile file of", nseq, "blocks into chunks of", chunksize, "blocks")
for i, batch in enumerate(batch_iterator(to_keep, chunksize)): # change
    filename = str(maf_folder) + "chr" + str(chrom) + "_%i.maf" % (i + 1)
    with open(filename, "w") as maf_handle:
        count = Bio.AlignIO.write(batch, maf_handle, "maf")
    print("Wrote %i sequences to %s" % (count, filename))

    filename = str(fasta_folder) + "chr" + str(chrom) + "_%i.fasta" % (i + 1)
    with open(filename, "w") as fasta_handle:
        count = Bio.AlignIO.write(batch, fasta_handle, "fasta")
    print("Wrote %i sequences to %s" % (count, filename))
