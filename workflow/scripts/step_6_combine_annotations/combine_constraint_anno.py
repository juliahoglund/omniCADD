import argparse
import glob
import os
import pandas as pd
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Combine constraint annotations for a chromosome.")
    parser.add_argument("-c", "--chromosome", type=str, required=True, help="chromosome number")
    parser.add_argument("-f", "--phast_folder", type=str, default="results/annotation/phast/phastCons", help="path to phastCons scores")
    parser.add_argument("-g", "--phylo_folder", type=str, default="results/annotation/phast/phyloP", help="path to phyloP scores")
    parser.add_argument("-i", "--gerp_folder", type=str, default="results/annotation/gerp", help="path to gerp scores")
    parser.add_argument("-o", "--output", type=str, default=None, help="output file name (default: constraint.<chromosome>.bed)")
    return parser.parse_args()

def load_and_concat(folder, chrom, suffix):
    # Try both: files directly in folder and in folder/chrN/
    files = []
    # 1. Check in folder/chrN/
    chr_folder = os.path.join(folder, f"chr{chrom}")
    if os.path.isdir(chr_folder):
        all_chr_files = os.listdir(chr_folder)
        if suffix == "phast.bed":
            files += [os.path.join(chr_folder, f) for f in all_chr_files if f.startswith(f"chr{chrom}_multiway-") and f.endswith(".phast.bed")]
        elif suffix == "phylo.bed":
            files += [os.path.join(chr_folder, f) for f in all_chr_files if f.startswith(f"chr{chrom}_multiway-") and f.endswith(".phylo.bed")]
    # 2. Also check in the main folder
    try:
        all_files = os.listdir(folder)
        if suffix == "phast.bed":
            files += [os.path.join(folder, f) for f in all_files if f.startswith(f"chr{chrom}_multiway-") and f.endswith(".phast.bed")]
        elif suffix == "phylo.bed":
            files += [os.path.join(folder, f) for f in all_files if f.startswith(f"chr{chrom}_multiway-") and f.endswith(".phylo.bed")]
    except FileNotFoundError:
        print(f"Directory {folder} does not exist.")
    # 3. Fallback for other suffixes (e.g. GERP)
    if not files and suffix not in ("phast.bed", "phylo.bed"):
        pattern = os.path.join(folder, f"chr{chrom}", f"chr{chrom}-*.{suffix}")
        print(f"Looking for files with pattern: {pattern}")
        files = sorted(glob.glob(pattern))
    # Sort files by chunk number (assumes format chrN_multiway-<chunk>.<suffix>)
    def extract_chunk_number(filename):
        import re
        match = re.search(r'_multiway-(\d+)\.', filename)
        return int(match.group(1)) if match else -1
    files = sorted(files, key=extract_chunk_number)
    print(f"Filtered and sorted files for {suffix}")
    if not files:
        print(f"No files found for {suffix} in {folder} or {chr_folder}")
    dfs = []
    for f in files:
        print(f"Reading file: {f}")
        df = pd.read_csv(f, sep="\t", header=None)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

def main():
    args = parse_args()
    chrom = args.chromosome

    print("Loading and combining phastCons files...")
    phastC = load_and_concat(args.phast_folder, chrom, "phast.bed")
    if phastC.empty:
        raise RuntimeError("No phastCons files found.")
    if phastC.shape[1] == 5:
        phastC = phastC.drop([0,3], axis=1)
        phastC.columns = ['start', 'end', 'phastCons']
    elif phastC.shape[1] == 3:
        phastC.columns = ['start', 'end', 'phastCons']
    elif phastC.shape[1] == 4:
        phastC = phastC.iloc[:, [1,2,3]]
        phastC.columns = ['start', 'end', 'phastCons']
    else:
        raise ValueError(f"Unexpected number of columns in phastCons input: {phastC.shape[1]}")

    print("Loading and combining phyloP files...")
    phyloP = load_and_concat(args.phylo_folder, chrom, "phylo.bed")
    if phyloP.empty:
        raise RuntimeError("No phyloP files found.")
    if phyloP.shape[1] == 5:
        phyloP = phyloP.drop([0,3], axis=1)
        phyloP.columns = ['start', 'end', 'phyloP']
    elif phyloP.shape[1] == 3:
        phyloP.columns = ['start', 'end', 'phyloP']
    elif phyloP.shape[1] == 4:
        phyloP = phyloP.iloc[:, [1,2,3]]
        phyloP.columns = ['start', 'end', 'phyloP']
    else:
        raise ValueError(f"Unexpected number of columns in phyloP input: {phyloP.shape[1]}")

    print("Merging phastCons and phyloP...")
    phast = pd.merge(phastC, phyloP, on=['start', 'end'], how='outer')
    phast['chr'] = f"chr{chrom}"
    phast = phast[['chr', 'start', 'end', 'phastCons', 'phyloP']]
    phast = phast.reset_index(drop=True)
    phast['Pos'] = phast.index + 1

    print("Finding and loading GERP file...")
    gerp_file = None
    for f in os.listdir(args.gerp_folder):
        if f"chr{chrom}" in f and f.endswith(".parsed"):
            gerp_file = os.path.join(args.gerp_folder, f)
            break
    if not gerp_file:
        raise RuntimeError(f"No GERP .parsed file found for chromosome {chrom} in {args.gerp_folder}")

    gerp = pd.read_csv(gerp_file, sep="\t", header=None, names=['Pos', 'GERP_NS', 'GERP_RS'])

    print("Merging all annotations on Pos...")
    constraint = pd.merge(phast, gerp, on='Pos', how='inner')
    # Drop start and end columns after merge
    constraint = constraint.drop(columns=['start', 'end'])
    # Rename chromosome column to #Chrom, remove 'chr' prefix, and make integer
    constraint = constraint.rename(columns={'chr': '#Chrom'})
    constraint['#Chrom'] = constraint['#Chrom'].astype(str).str.replace('chr', '', regex=False).astype(int)
    constraint = constraint[['#Chrom', 'Pos', 'GERP_NS', 'GERP_RS', 'phastCons', 'phyloP']]

    print("Writing output...")
    output_file = args.output or f"constraint.{chrom}.bed"
    constraint.to_csv(output_file, sep="\t", index=False, header=True, float_format="%.6f")
    print("Done.")

if __name__ == "__main__":
    main()
