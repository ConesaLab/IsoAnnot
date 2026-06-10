#!/usr/bin/env python3
"""
Author: Raúl López

Script to expand Nuclear Localization Signal (NLS) annotations from unique representative 
sequences back to all original isoforms. This script reverses the deduplication process 
by using a mapping file to ensure every isoform receives the functional annotation 
derived from its protein sequence.
"""

import argparse
import pandas as pd
import sys
import logging

def main():
    """
    Main execution logic to map NLS results back to original transcript identifiers.

    Args:
        --parsed (str): Path to the TSV file containing NLS results for unique IDs.
        --mapping (str): Path to the TSV mapping file (UniqueID vs OriginalIDs).
        --output (str): Path to save the final expanded TSV output.
    """
    parser = argparse.ArgumentParser(description="Expand NLS annotations from representatives to all isoforms")
    parser.add_argument("--parsed", required=True, help="Output from parse_nls.py (TSV)")
    parser.add_argument("--mapping", required=True, help="Mapping file from nls_deduplicate.py")
    parser.add_argument("--output", required=True, help="Final expanded NLS output (TSV)")
    args = parser.parse_args()

    # 1. Load the mapping file (Representative ID -> List of original IDs)
    logging.info(f"Loading mapping file from {args.mapping}")
    try:
        mapping_df = pd.read_csv(args.mapping, sep='\t')
        # Split the comma-separated string of original IDs into a Python list
        mapping_df['OriginalIDs'] = mapping_df['OriginalIDs'].str.split(',')
    except Exception as e:
        logging.error(f"Error reading mapping file: {e}")
        sys.exit(1)

    # 2. Convert to "long" format: one row per OriginalID
    # The 'explode' function creates a new row for each element in the OriginalIDs list,
    # effectively duplicating the UniqueID for each associated original transcript.
    logging.info("Exploding mapping to original isoform resolution")
    mapping_long = mapping_df.explode('OriginalIDs').rename(
        columns={'OriginalIDs': 'ID_original', 'UniqueID': 'ID'}
    )

    # 3. Load the NLS results (contains data for representative IDs only)
    logging.info(f"Reading NLS parsed results from {args.parsed}")
    try:
        nls_results = pd.read_csv(args.parsed, sep='\t')
    except Exception as e:
        logging.error(f"Error reading parsed NLS file: {e}")
        sys.exit(1)

    # 4. Merge: Combine NLS results with the expanded mapping
    # This step assigns the specific NLS metadata to every original protein ID 
    # that shared the same sequence during deduplication.
    final_df = pd.merge(mapping_long, nls_results, on='ID', how='inner')

    # 5. Final Cleanup: Standardize identifiers
    # Replace the representative 'ID' column with the actual original isoform ID
    final_df['ID'] = final_df['ID_original']
    # Remove the auxiliary column used for the join operation
    final_df = final_df.drop(columns=['ID_original'])

    # 6. Save the expanded results to a TSV file
    logging.info(f"Saving expanded annotations to {args.output}")
    try:
        final_df.to_csv(args.output, sep='\t', index=False, float_format='%.3f')
    except IOError as e:
        logging.error(f"Failed to write output file: {e}")
        sys.exit(1)

    logging.info(f"Expansion complete. Total annotations generated: {len(final_df)}")

if __name__ == "__main__":
    # Initialize logging with timestamp-free clean formatting
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
