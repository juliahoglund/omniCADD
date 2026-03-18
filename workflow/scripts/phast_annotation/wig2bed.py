import os
import pandas as pd
import argparse
from natsort import natsorted


def main():
    parser = argparse.ArgumentParser(description="Convert wig files to bed format.")
    parser.add_argument("-d", "--directory", type=str, default="results/annotation/phast/phastCons/chr1/",
                        help="Directory with phastCons / phyloP wig scores")
    parser.add_argument("-t", "--type", type=str, default="phast",
                        help="Type of analysis [phast | phylo]")
    args = parser.parse_args()

    files = natsorted([f for f in os.listdir(args.directory) if f.endswith(".wig")])
    last_row = 0

    for i, file in enumerate(files):
        chromosome = file.split('_')[0].replace('chr', '')
        part = pd.read_csv(os.path.join(args.directory, file), header=None, sep="\t")

        # Skip rows that contain header
        part = part[~part.iloc[:, 0].astype(str).str.startswith("fixedStep")].reset_index(drop=True)

        part['chromosome'] = chromosome

        # Update start positions to ensure continuity across files and make them 1-based
        part['start'] = range(last_row + 1, last_row + 1 + len(part))
        part['end'] = part['start'] + 1

        part['score'] = part.iloc[:, 0]
        part = part[['chromosome', 'start', 'end', 'score']]

        output_file = os.path.join(args.directory, f"{file.split('.')[0]}.{args.type}.bed")
        part.to_csv(output_file, header=False, index=False, sep="\t")

        # Update last_row to the last end position
        last_row = part['end'].iloc[-1]


if __name__ == "__main__":
    main()
