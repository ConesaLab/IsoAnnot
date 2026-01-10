#!/usr/bin/env python3
import argparse, csv, sys, logging
from IsoAnnot import openfile

def get_structural_classification_prot(classification_filename):
    """
    Reads the SQANTI classification file to create a Protein-to-Isoform map.
    Handles different SQANTI formats by checking for the existence of 'refProt'.
    """
    structural_classification = {}
    with openfile(classification_filename, "rt") as read_handler:
        tsv_reader = csv.DictReader(read_handler, delimiter="\t")
        for line in tsv_reader:
            # Key priority: refProt (if exists and is not NA), otherwise isoform ID
            ref_prot = line.get("refProt", None)
            iso_id = line.get("isoform")
            
            # Map based on the protein ID used by NucImport
            key = ref_prot if ref_prot and ref_prot != "NA" else iso_id
            structural_classification[key] = {"isoform": iso_id}
            
    return structural_classification

def main():
    """
    Maps consensus NLS predictions back to transcript/isoform IDs.
    """
    parser = argparse.ArgumentParser(description='IsoAnnot NLS Structural Layer')
    parser.add_argument('--consensus_file', required=True, help="Output from parse_nls.py")
    parser.add_argument('--classification_file', required=True, help="sqanti_classification.txt")
    parser.add_argument('--output', required=True, help="Output GTF file")
    args = parser.parse_args()

    # Load ID mapping from SQANTI
    transcript_info = get_structural_classification_prot(args.classification_file)

    with open(args.consensus_file, "r") as f_in, open(args.output, "w") as f_out:
        logging.info(f"Starting NLS conversion to GTF {args.output}")
        reader = csv.DictReader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")
        
        for row in reader:
            prot_id = row['ID']
            if prot_id in transcript_info:
                iso_id = transcript_info[prot_id]['isoform']
                
                attr = (f"ID=NLS_{row['Class']}; Name=NLS_{row['Class']}; Desc=NLS_NucImport")
                
                # Write GTF row: seqid, source, feature, start, end, score, strand, frame, attributes
                writer.writerow([
                    iso_id, 
                    "NucImport", 
                    "NLS",
                    row['Pos'], 
                    int(row['Pos']) + 20, # Fixed 20aa window
                    ".", ".", ".", attr
                ])

if __name__ == "__main__":
    main()
