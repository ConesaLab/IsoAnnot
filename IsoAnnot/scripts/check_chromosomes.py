#!/usr/bin/env python3
"""
Author: Raúl López

Script to verify chromosome naming consistency between a query transcriptome (GTF),
a reference GTF, and a reference FASTA file.
"""

import sys
import argparse
import gzip

def open_file(path):
    """
    Opens a file handles both gzipped and plain text files.

    Args:
        path (str): Path to the input file.

    Returns:
        file_object: A file handle in text read mode.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def get_ids_gtf(path):
    """
    Extracts unique chromosome IDs from the first column of a GTF file.

    Args:
        path (str): Path to the GTF file.

    Returns:
        set: A set containing all unique chromosome/scaffold names.
    """
    ids = set()
    with open_file(path) as f:
        for line in f:
            # Skip comments and empty lines
            if line.startswith("#") or not line.strip(): 
                continue
            # Chromosome ID is always the first column
            ids.add(line.split("\t")[0])
    return ids

def get_ids_fasta(path):
    """
    Extracts unique sequence IDs from a FASTA file header.

    Args:
        path (str): Path to the FASTA file.

    Returns:
        set: A set containing all sequence names (stripped of '>' and descriptions).
    """
    ids = set()
    with open_file(path) as f:
        for line in f:
            if line.startswith(">"):
                # Extract only the ID before the first whitespace
                ids.add(line[1:].split()[0])
    return ids

def main():
    """
    Main logic to compare chromosome ID sets and report inconsistencies or warnings.
    """
    try:
        parser = argparse.ArgumentParser(description='Chromosome Consistency Checker')
        parser.add_argument("--user_input", required=True, help="Query transcriptome (GTF or FASTA)")
        parser.add_argument("--ref_gtf", required=True, help="Reference GTF file")
        parser.add_argument("--ref_fasta", required=True, help="Reference FASTA genome")
        args = parser.parse_args()

        # Step 1: Handle input type. If user input is FASTA, bypass consistency check 
        # as coordinates are not relevant for naming match in the same way as GTF.
        if str(args.user_input).lower().endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz")):
            sys.exit(0)
        else:
            user_ids = get_ids_gtf(args.user_input)

        # Step 2: Load reference IDs
        ref_gtf_ids = get_ids_gtf(args.ref_gtf)
        ref_fasta_ids = get_ids_fasta(args.ref_fasta)

        # Step 3: Verify FASTA consistency (Physical presence)
        # All chromosomes in the user GTF MUST exist in the Reference FASTA
        missing_fasta = user_ids - ref_fasta_ids

        if missing_fasta:
            example_user = sorted(list(missing_fasta))
            example_ref = sorted(list(ref_fasta_ids))
            sys.stderr.write("CHROMOSOME CONSISTENCY ERROR: ")
            sys.stderr.write(f"[-] Transcriptome IDs: {example_user[:5]} ... NOT found in Reference FASTA: {example_ref[:5]} ...\n")
            sys.stderr.write("ADVICE: Ensure Transcriptome file and FASTA reference use the same convention (e.g., '1' vs 'chr1').\n")
            sys.exit(1)

        # Step 4: Verify GTF consistency (Annotation overlap)
        matches_gtf = user_ids & ref_gtf_ids
        missing_gtf = user_ids - ref_gtf_ids
        match_percentage = (len(matches_gtf) / len(user_ids)) * 100 if user_ids else 0

        # Fatal error if no chromosomes match the reference GTF
        if match_percentage == 0:
            example_missing = sorted(list(missing_gtf))
            example_ref_gtf = sorted(list(ref_gtf_ids))
            sys.stderr.write("CHROMOSOME CONSISTENCY ERROR: ")
            sys.stderr.write(f"Your IDs do not match the reference GTF.\n")
            sys.stderr.write(f"[-] Transcriptome IDs: {example_missing[:5]} NOT found in GTF: {example_ref_gtf[:5]}\n")
            sys.exit(1)

        # Warning if match is partial (typical when using primary assemblies with many scaffolds)
        elif match_percentage < 99:
            example_missing = sorted(list(missing_gtf))
            sys.stderr.write(f"\n[!] WARNING: Only {match_percentage:.2f}% of your transcriptome IDs are in the reference GTF.\n")
            sys.stderr.write(f"    Transcripts on unannotated scaffolds will be labeled as 'Novel' by SQANTI3.\n")
            sys.stderr.write(f"    Example of unannotated ID: {example_missing[0]}\n\n")
        
        sys.exit(0)

    except Exception as e:
        sys.stderr.write(f"FATAL ERROR: {str(e)}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
