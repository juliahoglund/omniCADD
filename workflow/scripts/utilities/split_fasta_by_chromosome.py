#!/usr/bin/env python3
"""
Split multi-FASTA file by chromosome

Takes a consensus FASTA with all chromosomes and splits into separate files.
Used in alignment generation workflow for GERP preprocessing.

Author: Julia Höglund
Date: 2026-02-27
"""

import sys
import os
from pathlib import Path
from Bio import SeqIO
from data_helper import split_fasta_by_chromosome

if __name__ == "__main__":
    # Get parameters from snakemake
    input_fasta = snakemake.input[0]
    output_dir = snakemake.params.outdir
    chromosomes = snakemake.params.chromosomes

    # Split FASTA
    split_fasta_by_chromosome(input_fasta, output_dir, chromosomes)
