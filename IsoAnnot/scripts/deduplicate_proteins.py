#!/usr/bin/env python3
"""
Author: Raúl López

Script to deduplicate protein sequences before NLS analysis. 
Identical sequences are grouped together, keeping one representative for processing 
to reduce computational cost while maintaining a mapping to all original IDs.
"""

import argparse
import sys
import logging
from Bio import SeqIO
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def main():
    """
    Main execution logic for sequence deduplication and mapping.
    
    Args:
        --input (str): Path to the input FASTA file.
        --output_fasta (str): Path to save unique representative sequences.
        --output_mapping (str): Path to save the TSV mapping (UniqueID vs OriginalIDs).
    """
    parser = argparse.ArgumentParser(description="Deduplicate protein sequences before NLS analysis")
    parser.add_argument("--input", required=True, help="Input fasta from nls_filter.py")
    parser.add_argument("--output_fasta", required=True, help="Fasta with unique representative sequences")
    parser.add_argument("--output_mapping", required=True, help="TSV mapping: UniqueID <tab> List of OriginalIDs")
    args = parser.parse_args()

    # Dictionary to store groups of IDs sharing the same sequence
    # Key: Sequence string, Value: List of protein identifiers
    seq_groups = defaultdict(list)

    # 1. Group IDs by identical sequence
    logging.info(f"Reading and grouping sequences from {args.input}")
    try:
        with open(args.input, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_groups[str(record.seq)].append(record.id)
    except Exception as e:
        logging.error(f"Failed to parse FASTA file: {e}")
        sys.exit(1)

    unique_records = []

    # 2. Create mapping file and select representative sequences
    logging.info(f"Generating mapping and unique records")
    try:
        with open(args.output_mapping, "w") as map_file:
            # Writing header for the mapping file
            map_file.write("UniqueID\tOriginalIDs\n")

            for seq, original_ids in seq_groups.items():
                # The first ID encountered is designated as the 'UniqueID' (Representative)
                representative_id = original_ids[0]

                # Save the relationship: Representative -> All original identifiers (comma separated)
                map_file.write(f"{representative_id}\t{','.join(original_ids)}\n")

                # Construct a new SeqRecord for the unique FASTA output
                unique_records.append(SeqRecord(Seq(seq), id=representative_id, description=""))
    except IOError as e:
        logging.error(f"Failed to write mapping file: {e}")
        sys.exit(1)

    # 3. Save FASTA containing only unique sequences
    try:
        with open(args.output_fasta, "w") as output_handle:
            SeqIO.write(unique_records, output_handle, "fasta")
    except IOError as e:
        logging.error(f"Failed to write unique FASTA file: {e}")
        sys.exit(1)

    # Summary statistics
    original_count = sum(len(v) for v in seq_groups.values())
    unique_count = len(unique_records)
    logging.info("Deduplication complete.")
    logging.info(f"Original proteins processed: {original_count}")
    logging.info(f"Unique representative proteins: {unique_count}")
    logging.info(f"Redundancy reduced by: {((original_count - unique_count) / original_count * 100):.2f}%")

if __name__ == "__main__":
    # Configure logging to output to stderr with clean formatting
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
