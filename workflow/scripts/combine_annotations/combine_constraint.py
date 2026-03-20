#!/usr/bin/env python3
"""Aggregate per-chunk GERP and PHAST scores into a chromosome-level constraint BED.

Reads per-chunk GERP rates.parsed, phastCons BED, phyloP BED and alignment
index files, then merges them into a single space-separated BED file suitable
for merge_annotations.py (-b / --bed argument).

Output columns: chr start end GERP_NS GERP_RS phastCons phyloP
"""

import argparse
import re
import sys
import pandas as pd


def read_index(path):
    """Parse alignment index file returning list of 0-based genomic start positions."""
    positions = []
    with open(path) as fh:
        for line in fh:
            m = re.match(r'start:\s*(\d+),\s*size:\s*(\d+)', line.strip())
            if m:
                start = int(m.group(1))
                size = int(m.group(2))
                positions.extend(range(start, start + size))
    return positions


def read_gerp(path, positions):
    """Read rates.parsed (GERP_NS\\tGERP_RS per position) aligned to genomic positions."""
    df = pd.read_csv(path, sep='\t', header=None, names=['GERP_NS', 'GERP_RS'])
    n = min(len(df), len(positions))
    df = df.iloc[:n].copy()
    df['start'] = positions[:n]
    return df[['start', 'GERP_NS', 'GERP_RS']]


def read_bed_score(path, col_name):
    """Read wig2bed BED file (0-based start), returning position + score columns."""
    try:
        df = pd.read_csv(path, sep='\t', header=None,
                         names=['chr', 'start', 'end', 'name', 'score'])
        return df.rename(columns={'score': col_name})[['start', col_name]]
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame(columns=['start', col_name])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gerp', nargs='+', required=True,
                        help='GERP rates.parsed files, one per chunk')
    parser.add_argument('--phastCons', nargs='+', required=True,
                        help='phastCons BED files (wig2bed output), one per chunk')
    parser.add_argument('--phyloP', nargs='+', required=True,
                        help='phyloP BED files (wig2bed output), one per chunk')
    parser.add_argument('--index', nargs='+', required=True,
                        help='Alignment index files, one per chunk')
    parser.add_argument('--chr', required=True, help='Chromosome identifier')
    parser.add_argument('-o', '--output', required=True,
                        help='Output constraint BED file (space-separated)')
    args = parser.parse_args()

    gerp_frames, phast_frames, phylo_frames = [], [], []

    for gerp_f, phast_f, phylo_f, idx_f in zip(
        sorted(args.gerp), sorted(args.phastCons),
        sorted(args.phyloP), sorted(args.index)
    ):
        positions = read_index(idx_f)
        if not positions:
            print(f"Warning: empty index file {idx_f}, skipping chunk", file=sys.stderr)
            continue
        gerp_frames.append(read_gerp(gerp_f, positions))
        phast_frames.append(read_bed_score(phast_f, 'phastCons'))
        phylo_frames.append(read_bed_score(phylo_f, 'phyloP'))

    if not gerp_frames:
        print("Error: no valid chunks found", file=sys.stderr)
        sys.exit(1)

    gerp_all = pd.concat(gerp_frames, ignore_index=True)
    phast_all = pd.concat(phast_frames, ignore_index=True) if phast_frames else pd.DataFrame(
        columns=['start', 'phastCons'])
    phylo_all = pd.concat(phylo_frames, ignore_index=True) if phylo_frames else pd.DataFrame(
        columns=['start', 'phyloP'])

    merged = gerp_all.merge(phast_all, on='start', how='outer')
    merged = merged.merge(phylo_all, on='start', how='outer')
    merged = merged.fillna(0)

    merged.insert(0, 'chr', f'chr{args.chr}')
    merged['end'] = merged['start']

    out_cols = ['chr', 'start', 'end', 'GERP_NS', 'GERP_RS', 'phastCons', 'phyloP']
    merged = merged[out_cols].sort_values('start')
    merged.to_csv(args.output, sep=' ', index=False)
    print(f"Wrote {len(merged)} positions to {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
