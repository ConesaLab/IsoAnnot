#!/usr/bin/env python3
"""
Author: Raúl López

Script to parse and filter predictions from NucImport output. 
It categorizes Nuclear Localization Signals (NLS) based on:
1. Import Probability (ImpProb): Probability of the protein being imported to the nucleus.
2. cNLS Score (cNLSProb): Probability of the presence of a classical NLS motif.

High-confidence predictions (meeting both thresholds) are saved to the output, 
while uncertain cases (high import but low motif score) are logged for manual review.
"""

import pandas as pd
import argparse
import logging
import sys

def main():
    """
    Main execution logic for parsing NucImport text results into a structured TSV.

    Args:
        --input (str): Path to the raw NucImport output file.
        --threshold_imp (float): Minimum probability for nuclear import (Default: 0.7).
        --threshold_cnls (float): Minimum probability for cNLS motif presence (Default: 0.3).
        --output (str): Path to save the filtered high-confidence NLS predictions.
    """
    parser = argparse.ArgumentParser(description="NLS Parser for NucImport results filtering")
    parser.add_argument("--input", required=True, help="NucImport raw output file")
    parser.add_argument("--threshold_imp", type=float, default=0.7, help="Minimum Import Probability threshold")
    parser.add_argument("--threshold_cnls", type=float, default=0.3, help="Minimum cNLS Probability threshold")
    parser.add_argument("--output", required=True, help="Path for High confidence NLS output (TSV)")
    args = parser.parse_args()

    # 1. Load the NucImport output
    # NucImport results usually have a non-standard whitespace separator and header stars
    logging.info(f"Reading NucImport results from {args.input}")
    try:
        # Use regex separator '\s+' to handle variable spaces between columns
        # Skip the first 2 lines (usually decorative star headers in NucImport)
        df = pd.read_csv(args.input, sep=r'\s+', skiprows=2, header=None,
                         names=['ID', 'ImpProb', 'cNLSProb', 'Seq', 'Class', 'Pos'])
    except Exception as e:
        logging.error(f"Error reading file {args.input}: {e}")
        sys.exit(1)

    if df.empty:
        logging.warning("The input file appears to be empty. No predictions to process.")
        sys.exit(1)

    # 2. Apply Filters based on biological confidence
    # High Confidence: Both Import and Motif probabilities meet the user thresholds
    mask_predicted = (df['ImpProb'] >= args.threshold_imp) & (df['cNLSProb'] >= args.threshold_cnls)
    df_nls_predicted = df[mask_predicted].copy()

    # Uncertain cases: High probability of being imported, but no clear cNLS motif found
    # These are worth logging as they might use non-classical import pathways
    mask_uncertain = (df['ImpProb'] >= args.threshold_imp) & (df['cNLSProb'] < args.threshold_cnls)
    df_uncertain = df[mask_uncertain].copy()

    # 3. Format and Save Output
    # Convert 'Pos' to integer to remove trailing decimals from the output file
    logging.info(f"Writing validated results to {args.output}")
    df_nls_predicted['Pos'] = df_nls_predicted['Pos'].astype(int)
    df_nls_predicted.to_csv(args.output, sep='\t', index=False, float_format='%.3f')

    # 4. Log Uncertain Cases for traceability
    # This provides the researcher with a list of proteins that might be nuclear 
    # despite lacking a strong cNLS signal.
    if not df_uncertain.empty:
        logging.warning("### UNCERTAIN NLS LOG (High Import / Low Motif Score) ###")
        logging.warning("ID\tImpProb\tcNLSProb\tPos\tClass")
        for _, row in df_uncertain.iterrows():
            logging.warning(f"{row['ID']}\t{row['ImpProb']:.3f}\t{row['cNLSProb']:.3f}\t{int(row['Pos'])}\t{row['Class']}")

    # Final summary for Snakemake logs
    logging.info(f"Parsing complete.")
    logging.info(f"High-confidence NLS (exported): {len(df_nls_predicted)}")
    logging.info(f"Uncertain NLS (logged as warnings): {len(df_uncertain)}")

if __name__ == "__main__":
    # Initialize logging to standard error for pipeline compatibility
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
