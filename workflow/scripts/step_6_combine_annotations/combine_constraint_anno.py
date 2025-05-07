#!/usr/bin/env python3

import os
import pandas as pd
import argparse
import glob
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Argument parser
parser = argparse.ArgumentParser(description="Combine constraint annotation files by chromosome.")
parser.add_argument("-c", "--chromosome", type=str, default="1", help="Chromosome number")
parser.add_argument("-f", "--phast-folder", type=str, default="results/annotation/phast/phastCons", help="Path to phastCons scores")
parser.add_argument("-g", "--phylo-folder", type=str, default="results/annotation/phast/phyloP", help="Path to phyloP scores")
parser.add_argument("-i", "--gerp-folder", type=str, default="results/annotation/gerp", help="Path to GERP scores")
parser.add_argument("-o", "--output-folder", type=str, default="results/combined", help="Output folder for combined files")

args = parser.parse_args()

# Ensure output folder exists
os.makedirs(args.output_folder, exist_ok=True)

# Check if output file already exists
output_file = os.path.join(args.output_folder, f"constraint.chr{args.chromosome}.bed")
if os.path.exists(output_file):
    logging.warning(f"Output file already exists and will be overwritten: {output_file}")

# Initialize empty DataFrames for combining
phast_combined = pd.DataFrame()
phylo_combined = pd.DataFrame()
gerp_combined = pd.DataFrame()

# Process phastCons files
logging.info(f"Processing phastCons files for chromosome {args.chromosome}...")
phast_folder = os.path.join(args.phast_folder, f"chr{args.chromosome}")
if os.path.exists(phast_folder):
    phast_files = [os.path.join(phast_folder, f) for f in os.listdir(phast_folder) if f.endswith(".phast.bed")]
    if phast_files:
        try:
            phast_combined = pd.concat(
                [pd.read_csv(f, sep="\t", header=None, usecols=[1, 2, 4], names=["start", "end", "phastCons"]) for f in phast_files],
                ignore_index=True
            )
            phast_combined["start"] = phast_combined["start"].astype(int)
            phast_combined["end"] = phast_combined["end"].astype(int)
        except Exception as e:
            logging.error(f"Error processing phastCons files: {e}")
    else:
        logging.warning(f"No phastCons files found in folder: {phast_folder}")
else:
    logging.warning(f"PhastCons folder not found: {phast_folder}")

# Process phyloP files
logging.info(f"Processing phyloP files for chromosome {args.chromosome}...")
phylo_folder = os.path.join(args.phylo_folder, f"chr{args.chromosome}")
if os.path.exists(phylo_folder):
    phylo_files = [os.path.join(phylo_folder, f) for f in os.listdir(phylo_folder) if f.endswith(".phylo.bed")]
    if phylo_files:
        try:
            phylo_combined = pd.concat(
                [pd.read_csv(f, sep="\t", header=None, usecols=[1, 2, 4], names=["start", "end", "phyloP"]) for f in phylo_files],
                ignore_index=True
            )
            phylo_combined["start"] = phylo_combined["start"].astype(int)
            phylo_combined["end"] = phylo_combined["end"].astype(int)
        except Exception as e:
            logging.error(f"Error processing phyloP files: {e}")
    else:
        logging.warning(f"No phyloP files found in folder: {phylo_folder}")
else:
    logging.warning(f"PhyloP folder not found: {phylo_folder}")

# Process GERP file
logging.info(f"Processing GERP file for chromosome {args.chromosome}...")
gerp_pattern = f"*chr{args.chromosome}_multiway.fa.rates"
gerp_files = glob.glob(os.path.join(args.gerp_folder, gerp_pattern))
if gerp_files:
    try:
        gerp_file = gerp_files[0]
        gerp_combined = pd.read_csv(gerp_file, sep="\t", header=None, names=["GERP_NS", "GERP_RS"])
        gerp_combined["start"] = range(1, len(gerp_combined) + 1)
        gerp_combined["chr"] = f"chr{args.chromosome}"
    except Exception as e:
        logging.error(f"Error processing GERP file: {e}")
else:
    logging.warning(f"GERP file not found for chromosome {args.chromosome}. Skipping GERP processing.")

# Combine phastCons and phyloP
logging.info("Combining phastCons and phyloP data...")
if not phast_combined.empty and not phylo_combined.empty:
    phast_phylo_combined = pd.merge(phast_combined, phylo_combined, on=["start", "end"], how="outer")
    phast_phylo_combined["chr"] = f"chr{args.chromosome}"
else:
    phast_phylo_combined = pd.DataFrame(columns=["chr", "start", "end", "phyloP", "phastCons"])

# Combine all data
logging.info("Combining all data...")
if not gerp_combined.empty and not phast_phylo_combined.empty:
    constraint_combined = pd.merge(gerp_combined, phast_phylo_combined, on=["chr", "start"], how="outer")
else:
    constraint_combined = pd.DataFrame(columns=["chr", "start", "end", "GERP_NS", "GERP_RS", "phyloP", "phastCons"])

# Reorder columns to: chr, start, end, GERP_NS, GERP_RS, phyloP, phastCons
constraint_combined = constraint_combined[["chr", "start", "end", "GERP_NS", "GERP_RS", "phyloP", "phastCons"]]

# Write the combined file
try:
    constraint_combined.to_csv(output_file, sep="\t", index=False, header=True)
    logging.info(f"Combined file written to {output_file}")
except Exception as e:
    logging.error(f"Error writing combined file: {e}")