#!/usr/bin/env python
# -*- coding: utf-8 -*-


# derived from:
# split_fasta.py (assumes you have biopython installed, e.g. with pip install biopython)
# from: https://thevirtuallaboratory.com/blog/splitting-a-multi-fasta
# Author: Bram van Dijk

# dependencies
import sys, math
import gzip
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
        sys.exit("usage: convert_alignments.py MAF_FILE N_CHUNKS OUTPUT_FOLDE_MAF REF_SPECIES")

mfile = sys.argv[1]  # maf file
ofile = gzip.open(mfile, "rt") \
                if mfile.endswith('.gz') else open(mfile, "r") # maf file


chunks=sys.argv[2] # number of chunks
maf_folder=sys.argv[3] # folder to save maf chunks in
ref_species=sys.argv[4]

to_keep = []
start = []
end = []

# parse alignment
for alignment in AlignIO.parse(ofile, "maf"):
    if ref_species in str(alignment):
        to_keep.append(alignment)
        for seqrec in alignment:
            if ref_species in seqrec.id:
                start.append(seqrec.annotations["start"])
                end.append(seqrec.annotations["start"]+seqrec.annotations["size"]-1) # next start will be this end if not -1
                # 0 based 1 based double check later


nseq = len(to_keep)
chunksize=math.ceil(nseq/int(chunks))
chrom = mfile.split('.')[0].split('chr')[1]

print("creatung indices for base pair positionvstart and end per maf block.")
for i, batch in enumerate(batch_iterator(start, chunksize)):
    indexfile = str(fasta_folder) + "/chr" + str(chrom) + "_%i.start" % (i +1)
    with open(indexfile, "w") as index_handle:
        index_handle.write('\n'.join([str(line) for line in batch]))

for i, batch in enumerate(batch_iterator(end, chunksize)):
    indexfile = str(fasta_folder) + "/chr" + str(chrom) + "_%i.stop" % (i +1)
    with open(indexfile, "w") as index_handle:
        index_handle.write('\n'.join([str(line) for line in batch]))

print("Splitting maffile file of", nseq, "blocks into chunks of", chunksize, "blocks")
for i, batch in enumerate(batch_iterator(to_keep, chunksize)):

    filename = str(maf_folder) + "/chr" + str(chrom) + "_%i.maf" % (i + 1)
    with open(filename, "w") as maf_handle:
        count = Bio.AlignIO.write(batch, maf_handle, "maf")
    print("Wrote %i sequences to %s" % (count, filename))
