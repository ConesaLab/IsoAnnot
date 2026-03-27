#!/usr/bin/env python3
import argparse
import pandas as pd
import sys
import logging

def main():
    parser = argparse.ArgumentParser(description="Expand NLS annotations from representatives to all isoforms")
    parser.add_argument("--parsed", required=True, help="Output from parse_nls.py (TSV)")
    parser.add_argument("--mapping", required=True, help="Mapping file from nls_deduplicate.py")
    parser.add_argument("--output", required=True, help="Final expanded NLS output (TSV)")
    args = parser.parse_args()

    # 1. Load the mapping file (Representative ID -> List of original IDs)
    try:        
        mapping_df = pd.read_csv(args.mapping, sep='\t')
        mapping_df['OriginalIDs'] = mapping_df['OriginalIDs'].str.split(',')
    except Exception as e:
        logging.error(f"Error reading mapping file: {e}")
        sys.exit(1)

    # 2. Convert to "long" format: one row per OriginalID
    # The 'explode' function creates a new row for each element in the OriginalIDs list
    mapping_long = mapping_df.explode('OriginalIDs').rename(columns={'OriginalIDs': 'ID_original', 'UniqueID': 'ID'})

    # 3. Load the NLS results (contains data for representative IDs only)
    nls_results = pd.read_csv(args.parsed, sep='\t')

    # 4. Merge: Combine NLS results with the expanded mapping
    # This assigns the NLS row to every original protein that shared that sequence
    final_df = pd.merge(mapping_long, nls_results, on='ID', how='inner')

    # 5. Final Cleanup: Replace the representative ID with the original isoform ID
    final_df['ID'] = final_df['ID_original']
    # Drop the auxiliary column used for mapping
    final_df = final_df.drop(columns=['ID_original'])

    # 6. Save the expanded results to a TSV file
    final_df.to_csv(args.output, sep='\t', index=False, float_format='%.3f')

    logging.info(f"Expansion complete. Total annotations: {len(final_df)}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
