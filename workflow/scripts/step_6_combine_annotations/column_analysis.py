#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Job van Schipstal
:Contact: job.vanschipstal@wur.nl
:Date: 18-02-2023
:Usage: see --help

Assess column correlation and relevance.
"""

# Import dependencies
from argparse import ArgumentParser
import pandas

# Create and process cli
parser = ArgumentParser(description=__doc__)
parser.add_argument("-d", "--derived",
                    help="annotated infile(s) containing derived variants to process", 
                    nargs="+", 
                    required=True)
parser.add_argument("-s", "--simulated",
                    help="annotated infile(s) containing simulated variants to process", 
                    nargs="+", 
                    required=True)
parser.add_argument("-o", "--outfolder",
                    help="Output folder for correlation tsv files",
                    type=str, 
                    default="None")

args = parser.parse_args()


def get_dataset(infiles: list):
    dataset = []
    for infile in infiles:
        dataset.append(pandas.read_csv(infile, sep='\t', na_values=['-'], low_memory=False))
    if len(dataset) > 1:
        return pandas.concat(dataset, axis=0)
    return dataset[0]


def analyse(dataset, name, out_folder):
    dataset.corr(numeric_only=True).to_csv(f"{out_folder}{name}_corr.tsv",
                                           sep="\t")
    return 1 - dataset.isnull().sum() / len(dataset)


out_folder = args.outfolder
if not out_folder.endswith("/"):
    out_folder += "/"

derived = get_dataset(args.derived)
derived.insert(0, "y", [0.0] * len(derived))
derived_rel = analyse(derived, "derived_variants", out_folder)

simulated = get_dataset(args.simulated)
simulated.insert(0, "y", [1.0] * len(simulated))
simulated_rel = analyse(simulated, "simulated_variants", out_folder)

combined = pandas.concat([simulated, derived], axis=0)
del simulated, derived
combined_rel = analyse(combined, "combined_variants", out_folder)

relevance_df = {"name": combined.columns}
with open(out_folder + "relevance.tsv", "w") as rel_file:
    rel_file.write("Column name\tDerived variants\t"
                   "Simulated variants\tCombined variants\n")
    for col, der_rel, sim_rel, com_rel in zip(combined.columns,
                                              derived_rel,
                                              simulated_rel,
                                              combined_rel):
        rel_file.write(f"{col}\t{der_rel:.3f}\t{sim_rel:.3f}\t{com_rel:.3f}\n")
