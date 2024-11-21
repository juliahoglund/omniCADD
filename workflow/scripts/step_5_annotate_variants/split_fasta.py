# dependencies
import os, sys, math
from argparse import ArgumentParser
from collections import defaultdict
from itertools import zip_longest

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def grouper(n, iterable):
    """ s -> (s0,s1,...sn-1), (sn,sn+1,...s2n-1), (s2n,s2n+1,...s3n-1), ... """
    FILLER = object()  # Value that couldn't be in data.
    for result in zip_longest(*[iter(iterable)]*n, fillvalue=FILLER):
        yield ''.join(v for v in result if v is not FILLER)

ffile = sys.argv[1] # fasta file to be formatted
chunks = sys.argv[2]

print("splitting", ffile, "into", chunks, "chunks")

species = defaultdict(list)
entries = defaultdict(list)

# retrieve all species present in any block to get maximum number of aligned species
for record in SeqIO.parse(ffile, "fasta"):
    record.id = record.id.split('.')[0]
    if record.id not in species:
        species[record.id] = []
    entries[record.id].append(len(record))
    species[record.id].append(record.seq)

ik = dict()
for i, k in enumerate(filedata):
    ik[i] = k   # dictionary key_of_index

nseq = entries[ik[0]][0]
chunksize=math.ceil(int(nseq)/int(chunks))
chrom = ffile.split('.')[0].split('chr')[1]

print("Splitting fasta file of", nseq, "basepairs into chunks of", chunksize, "bps per file")

filenames = []
for i in range(1, chunks+1):
    filenames.append("chr" + str(chrom) + "_multiway-%i.fa" % i)
filedata = {filename: open(filename, 'w') for filename in filenames}

for i, (key, value) in enumerate(species.items()):
    j = 0 
    for group in grouper(chunksize, str(value[0])):
        j += 1
        filedata[ik[j-1]].write('>' + key + '\n' + group + '\n')

for file in filedata.values():
    file.close()

