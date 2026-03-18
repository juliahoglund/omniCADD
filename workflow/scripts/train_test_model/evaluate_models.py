#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
:Author: Julia Höglund
:Date: 2024
:Usage: evaluate_models.py --help

Evaluates trained models by parsing their statistics files.
Compares different hyperparameter combinations (C and max_iter)
and identifies the best performing model based on ROC-AUC score.
Outputs a summary table and saves the best parameters to JSON.
"""

import sys
import json
import logging
import re
from argparse import ArgumentParser
from pathlib import Path
import pandas as pd

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = ArgumentParser(description=__doc__)
parser.add_argument("--stats",
                    help="Stats files from trained models (.stats.txt)",
                    type=str,
                    required=True,
                    nargs="+")

parser.add_argument("--summary",
                    help="Output file for summary table (TSV format)",
                    type=str,
                    required=True)

parser.add_argument("--best-params",
                    help="Output file for best parameters (JSON format)",
                    type=str,
                    required=True)


def parse_stats_file(stats_file):
    """
    Parse a stats file and extract key metrics.
    
    :param stats_file: Path to .stats.txt file
    :return: Dictionary with parsed metrics
    """
    try:
        with open(stats_file, 'r') as f:
            content = f.read()
        
        # Extract hyperparameters from filename or first line
        filename = Path(stats_file).stem
        
        # Try to extract from filename pattern: fold_X_YC_Ziter.mod.stats
        c_match = re.search(r'_(\d+\.?\d*)C_', filename)
        iter_match = re.search(r'_(\d+)iter', filename)
        fold_match = re.search(r'fold_(\d+)', filename)
        
        if not c_match or not iter_match:
            # Try to extract from file content
            header_match = re.search(r'Stats for c=(\d+\.?\d*), max_iter=(\d+)', content)
            if header_match:
                c = float(header_match.group(1))
                max_iter = int(header_match.group(2))
            else:
                logging.warning(f"Could not parse hyperparameters from {stats_file}")
                return None
        else:
            c = float(c_match.group(1))
            max_iter = int(iter_match.group(1))
        
        fold = int(fold_match.group(1)) if fold_match else None
        
        # Extract ROC-AUC score
        auc_match = re.search(r'ROC-AUC score:\s+(\d+\.?\d*)', content)
        if auc_match:
            auc = float(auc_match.group(1))
        else:
            logging.warning(f"Could not find ROC-AUC score in {stats_file}")
            auc = None
        
        # Extract confusion matrix values
        # Format: [[TN FP]
        #          [FN TP]]
        cm_match = re.search(r'\[\[(\d+)\s+(\d+)\]\s+\[(\d+)\s+(\d+)\]\]', content)
        if cm_match:
            tn, fp, fn, tp = map(int, cm_match.groups())
            accuracy = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) > 0 else None
            precision = tp / (tp + fp) if (tp + fp) > 0 else None
            recall = tp / (tp + fn) if (tp + fn) > 0 else None
        else:
            logging.warning(f"Could not parse confusion matrix from {stats_file}")
            accuracy = precision = recall = None
        
        return {
            'file': str(stats_file),
            'fold': fold,
            'c': c,
            'max_iter': max_iter,
            'auc': auc,
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall
        }
    
    except Exception as e:
        logging.error(f"Error parsing {stats_file}: {e}")
        return None


def main():
    args = parser.parse_args()
    
    try:
        logging.info(f"Parsing {len(args.stats)} stats files...")
        
        # Parse all stats files
        results = []
        for stats_file in args.stats:
            parsed = parse_stats_file(stats_file)
            if parsed:
                results.append(parsed)
        
        if not results:
            logging.error("No stats files could be parsed successfully")
            sys.exit(1)
        
        # Create DataFrame
        df = pd.DataFrame(results)
        
        # Calculate mean metrics across folds for each hyperparameter combination
        summary = df.groupby(['c', 'max_iter']).agg({
            'auc': ['mean', 'std'],
            'accuracy': ['mean', 'std'],
            'precision': ['mean', 'std'],
            'recall': ['mean', 'std']
        }).reset_index()
        
        # Flatten column names
        summary.columns = ['_'.join(col).strip('_') for col in summary.columns.values]
        
        # Sort by mean AUC (descending)
        summary = summary.sort_values('auc_mean', ascending=False)
        
        # Save summary table
        logging.info(f"Saving summary to {args.summary}")
        summary.to_csv(args.summary, sep='\t', index=False, float_format='%.6f')
        
        # Identify best parameters
        best_row = summary.iloc[0]
        best_params = {
            'c': float(best_row['c']),
            'max_iter': int(best_row['max_iter']),
            'mean_auc': float(best_row['auc_mean']),
            'std_auc': float(best_row['auc_std']) if pd.notna(best_row['auc_std']) else 0.0
        }
        
        # Save best parameters
        logging.info(f"Saving best parameters to {args.best_params}")
        logging.info(f"Best parameters: C={best_params['c']}, max_iter={best_params['max_iter']}, "
                     f"AUC={best_params['mean_auc']:.6f} ± {best_params['std_auc']:.6f}")
        
        with open(args.best_params, 'w') as f:
            json.dump(best_params, f, indent=2)
        
        logging.info("Model evaluation completed successfully")
        
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
