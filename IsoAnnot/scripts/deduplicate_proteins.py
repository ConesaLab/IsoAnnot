#!/usr/bin/env python3
import argparse
import sys
import logging
from Bio import SeqIO
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def main():
    parser = argparse.ArgumentParser(description="Deduplicate protein sequences before NLS analysis")
    parser.add_argument("--input", required=True, help="Input fasta from nls_filter.py")
    parser.add_argument("--output_fasta", required=True, help="Fasta with unique representative sequences")
    parser.add_argument("--output_mapping", required=True, help="TSV mapping: UniqueID <tab> List of OriginalIDs")
    args = parser.parse_args()

    seq_groups = defaultdict(list)

    # 1. Group IDs by identical sequence
    # Using a dictionary where key = sequence string and value = list of headers
    with open(args.input, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_groups[str(record.seq)].append(record.id)

    unique_records = []

    # 2. Create mapping file and select representatives
    with open(args.output_mapping, "w") as map_file:
        map_file.write("UniqueID\tOriginalIDs\n")

        for seq, original_ids in seq_groups.items():
            # Use the first ID in the list as the "UniqueID" (Representative)
            representative_id = original_ids[0]

            # Save relationship: Representative -> All identical sequences
            map_file.write(f"{representative_id}\t{','.join(original_ids)}\n")

            # Create the record for the new unique FASTA file
            unique_records.append(SeqRecord(Seq(seq), id=representative_id, description=""))

    # 3. Save FASTA of unique sequences
    with open(args.output_fasta, "w") as output_handle:
        SeqIO.write(unique_records, output_handle, "fasta")

    logging.info(f"Deduplication complete.")
    logging.info(f"Original proteins: {sum(len(v) for v in seq_groups.values())}")
    logging.info(f"Unique proteins to process: {len(unique_records)}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
