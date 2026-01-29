#!/usr/bin/env python3
import pandas as pd
import argparse
import logging
import sys

def main():
    """
    Processes a NucImport output file to filter NLS predictions.
    Filters results based on Import and cNLS probability thresholds.
    """
    parser = argparse.ArgumentParser(description="NLS Parser for NucImport")
    parser.add_argument("--input", required=True, help="NucImport output file")
    parser.add_argument("--threshold_imp", type=float, default=0.7, help="Minimum Import Probability")
    parser.add_argument("--threshold_cnls", type=float, default=0.3, help="Minimum cNLS Probability")
    parser.add_argument("--output", required=True, help="High confidence NLS output (TSV)")
    args = parser.parse_args()

    # 1. Load the output
    try:
        # NucImport format uses variable whitespace; skip first 2 star-header lines
        df = pd.read_csv(args.input, sep=r'\s+', skiprows=2, header=None,
                         names=['ID', 'ImpProb', 'cNLSProb', 'Seq', 'Class', 'Pos'])
    except Exception as e:
        logging.error(f"Error reading file {args.input}: {e}")
        sys.exit(1)

    if df.empty:
        logging.error("The input file is empty. Exiting.")
        sys.exit(1)

    # 2. Apply Filters
    # High Confidence: Both thresholds met
    mask_predicted = (df['ImpProb'] >= args.threshold_imp) & (df['cNLSProb'] >= args.threshold_cnls)
    df_nls_predicted = df[mask_predicted].copy()

    # Uncertain: High Import but low cNLS score
    mask_uncertain = (df['ImpProb'] >= args.threshold_imp) & (df['cNLSProb'] < args.threshold_cnls)
    df_uncertain = df[mask_uncertain].copy()

    # 3. Format and Save Output
    # Ensure Position is integer for clean output
    df_nls_predicted['Pos'] = df_nls_predicted['Pos'].astype(int)
    df_nls_predicted.to_csv(args.output, sep='\t', index=False, float_format='%.3f')

    # 4. Log Uncertain Cases to stderr/log
    if not df_uncertain.empty:
        logging.warning("### UNCERTAIN NLS LOG ###")
        logging.warning("Description: Proteins with high Import Probability but low cNLS motif score.")
        logging.warning("ID\tImpProb\tcNLSProb\tPos\tClass")
        for _, row in df_uncertain.iterrows():
            logging.warning(f"{row['ID']}\t{row['ImpProb']:.3f}\t{row['cNLSProb']:.3f}\t{int(row['Pos'])}\t{row['Class']}")

    logging.info(f"Processing finished.")
    logging.info(f"NLS predicted (saved to output): {len(df_nls_predicted)}")
    logging.info(f"Uncertain NLS (logged as warnings): {len(df_uncertain)}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
