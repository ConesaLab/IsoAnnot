#!/usr/bin/env python3
import pandas as pd
import argparse
import logging
import sys

def get_safe_mode(x):
    """
    Retrieves the most frequent value (mode) in a series safely.
    Returns the first mode if exists, otherwise "N/A".
    """
    m = x.mode()
    return m.iloc[0] if not m.empty else "N/A"

def main():
    """
    Processes 6 NucImport output files to generate a consensus annotation.
    Filters results based on Import and cNLS probability thresholds.
    """
    parser = argparse.ArgumentParser(description="Consensus NLS Parser for IsoAnnot")
    parser.add_argument("--inputs", nargs='+', required=True, help="List of 6 NucImport output files")
    parser.add_argument("--threshold_imp", type=float, default=0.7, help="Minimum average Import Probability")
    parser.add_argument("--threshold_cnls", type=float, default=0.5, help="Minimum average cNLS Probability")
    parser.add_argument("--output", required=True, help="High confidence NLS output (TSV)")
    args = parser.parse_args()

    all_dfs = []
    # Load all 6 model outputs
    for f in args.inputs:
        try:
            # NucImport format uses variable whitespace; skip first 2 star-header lines
            df = pd.read_csv(f, sep=r'\s+', skiprows=2, header=None,
                             names=['ID', 'ImpProb', 'cNLSProb', 'Seq', 'Class', 'Pos'])
            all_dfs.append(df)
        except Exception as e:
            logging.error(f"Error reading file {f}: {e}")
            continue
    
    if not all_dfs:
        logging.error("No valid data loaded. Exiting.")
        sys.exit(1)

    # Merge results from all models
    combined = pd.concat(all_dfs)

    # Aggregate: Mean for probabilities, Mode for categorical and positional data
    consensus = combined.groupby('ID').agg({
        'ImpProb': ['mean', 'std'],
        'cNLSProb': ['mean', 'std'],
        'Pos': get_safe_mode,
        'Seq': get_safe_mode,
        'Class': get_safe_mode
    })

    # Flatten multi-index columns
    consensus.columns = ['ImpProb_mean', 'ImpProb_std', 'cNLSProb_mean', 'cNLSProb_std', 'Pos', 'Seq', 'Class']
    consensus = consensus.reset_index()

    # Apply Confidence filter (High Import AND High cNLS)
    nls_predicted = (consensus['ImpProb_mean'] >= args.threshold_imp) & \
                (consensus['cNLSProb_mean'] >= args.threshold_cnls)
    
    df_nls_predicted = consensus[nls_predicted].copy()
    # Ensure Position is integer
    df_nls_predicted['Pos'] = df_nls_predicted['Pos'].astype(int)

    # Apply Uncertain filter (High Import BUT Weak cNLS)
    nls_uncertain = (consensus['ImpProb_mean'] >= args.threshold_imp) & \
                     (consensus['cNLSProb_mean'] < args.threshold_cnls)
    
    df_uncertain = consensus[nls_uncertain].copy()

    # Save outputs as TSV
    df_nls_predicted.to_csv(args.output, sep='\t', index=False, float_format='%.3f')
    logging.warning("### UNCERTAIN NLS LOG ###\n")
    logging.warning("Description: Proteins with high Import Probability but low cNLS motif score.\n")
    logging.warning("ID\tImpProb_mean\tImpProb_std\tcNLSProb_mean\tcNLSProb_std\tSuggested_Pos\tSuggested_Class\n")
    for _, row in df_uncertain.iterrows():
        logging.warning(f"{row['ID']}\t{row['ImpProb_mean']:.3f}\t{row['ImpProb_std']:.3f}\t"
                       f"{row['cNLSProb_mean']:.3f}\t{row['cNLSProb_std']:.3f}\t{row['Pos']}\t{row['Class']}\n")

    logging.info(f"Consensus completed. NLS predicted: {len(df_nls_predicted)}, Uncertain NLS: {len(df_uncertain)}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
