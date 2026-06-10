#!/usr/bin/env python3
"""
Author: Raúl López

Script to filter and pre-process protein sequences for Nuclear Localization Signal (NLS) prediction.
This module performs two critical tasks:
1. Length Filtering: Removes sequences shorter than 30 amino acids, as they are unlikely 
   to contain complex localization signals.
2. Sequence Cleaning: Replaces ambiguous amino acids (U, O, B, Z) with '@' to ensure 
   compatibility with downstream prediction tools like NucImport.
"""

import argparse, sys, os, re, logging, traceback
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def protein_dict(fasta_prot):
    """
    Parses a FASTA file, applies a length filter, and cleans ambiguous amino acids.

    Args:
        fasta_prot (str): Path to the input FASTA file containing protein sequences.

    Returns:
        dict: A dictionary where keys are protein IDs and values are cleaned sequence strings.
    """
    prots = dict()
    
    logging.info(f"Parsing and filtering sequences from: {fasta_prot}")
    
    # Use Biopython SeqIO to parse the FASTA file
    with open(fasta_prot, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            
            # 1. Length filter: skip proteins shorter than 30 aa
            # NLS signals usually require a minimum context to be predicted accurately
            if len(record.seq) < 30:
                continue

            # 2. Sequence cleaning:
            # Prediction tools often crash with non-standard amino acids.
            # We replace Selenocysteine (U), Pyrrolysine (O), and ambiguous 
            # asparagine/aspartic acid (B) or glutamine/glutamic acid (Z) with '@'.
            clean_seq = re.sub(r"[UOBZ]", "@", str(record.seq))

            # Store cleaned sequence in the dictionary
            prots[record.id] = clean_seq

    return prots

def main():
    """
    Main execution logic for the NLS protein pre-processing filter.
    """
    try:
        parser = argparse.ArgumentParser(description="NLS Filter and Pre-processing")
        parser.add_argument("--input", required=True, help="Input protein FASTA file")
        parser.add_argument("--output", required=True, help="Output filtered FASTA file")
        args = parser.parse_args()

        # Step 1: Process and filter the proteins
        prot_dict = protein_dict(args.input)

        # Step 2: Convert dictionary back to SeqRecord objects for writing
        filtered_records = []
        for p_id, p_seq in prot_dict.items():
            # We keep descriptions empty to ensure clean FASTA headers
            record = SeqRecord(Seq(p_seq), id=p_id, description="")
            filtered_records.append(record)

        # Step 3: Write the results to the output FASTA file
        logging.info(f"Writing {len(filtered_records)} filtered proteins to {args.output}")
        with open(args.output, "w") as output_handle:
            SeqIO.write(filtered_records, output_handle, "fasta")

        logging.info("Filtering and pre-processing complete.")

    except Exception as ex:
        logging.error(f"Error during NLS filtering: {str(ex)}")
        traceback.print_exc()
        # Return software error exit code
        sys.exit(os.EX_SOFTWARE)

if __name__ == "__main__":
    # Configure logging to provide clean feedback in Snakemake logs
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
