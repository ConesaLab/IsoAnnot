#!/usr/bin/env python3
"""
Author: Raúl López

Script to map Nuclear Localization Signal (NLS) predictions back to transcript IDs 
and format them into a GTF file. It utilizes the structural classification from 
SQANTI3 to link protein-level predictions with their parent isoforms.
"""

import argparse, csv, sys, logging
from IsoAnnot import openfile

def get_structural_classification_prot(classification_filename):
    """
    Parses the SQANTI classification file to build a lookup table between 
    protein identifiers and isoform IDs.

    Args:
        classification_filename (str): Path to the sqanti_classification.txt file.

    Returns:
        dict: A mapping where keys are the protein IDs (refProt or isoform) 
              and values are dictionaries containing the isoform ID.
    """
    structural_classification = {}
    
    with openfile(classification_filename, "rt") as read_handler:
        tsv_reader = csv.DictReader(read_handler, delimiter="\t")
        for line in tsv_reader:
            # Logic for ID matching:
            # 1. Check if 'refProt' exists (standard for reference-based analysis).
            # 2. Fallback to 'isoform' ID if refProt is missing or NA.
            ref_prot = line.get("refProt", None)
            iso_id = line.get("isoform")

            # Determine the key used by NucImport/downstream tools for this protein
            key = ref_prot if ref_prot and ref_prot != "NA" else iso_id
            structural_classification[key] = {"isoform": iso_id}

    return structural_classification

def main():
    """
    Main logic to translate NLS results into GTF format.
    Coordinates are kept at protein level (amino acids).
    """
    parser = argparse.ArgumentParser(description='IsoAnnot NLS Structural Layer Mapping')
    parser.add_argument('--nls_file', required=True, help="Parsed results from NucImport (TSV)")
    parser.add_argument('--classification_file', required=True, help="SQANTI classification file path")
    parser.add_argument('--output', required=True, help="Path for the output GTF file")
    args = parser.parse_args()

    # Configure basic logging for pipeline status
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Step 1: Load the protein-to-isoform mapping
    logging.info("Loading ID mapping from SQANTI classification...")
    transcript_info = get_structural_classification_prot(args.classification_file)

    # Step 2: Convert TSV predictions to GTF records
    with open(args.nls_file, "r") as f_in, open(args.output, "w") as f_out:
        logging.info(f"Starting NLS conversion to GTF: {args.output}")
        
        reader = csv.DictReader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")

        for row in reader:
            prot_id = row['ID']
            
            # Map protein prediction back to its parent isoform
            if prot_id in transcript_info:
                iso_id = transcript_info[prot_id]['isoform']

                # Format GTF attributes column
                attr = (f"ID=NLS_{row['Class']}; Name=NLS_{row['Class']}; Desc=NLS_NucImport")

                # Step 3: Write GTF row (9 columns standard)
                writer.writerow([
                    iso_id,              # seqid (Transcript ID)
                    "NucImport",         # source
                    "NLS",               # feature
                    row['Pos'],          # start (1-based AA position)
                    int(row['Pos']) + 20, # end (window specified by NucImport)
                    ".",                 # score
                    ".",                 # strand (not applicable in protein coords)
                    ".",                 # frame
                    attr                 # attributes
                ])

    logging.info("NLS layer generation completed.")

if __name__ == "__main__":
    main()
