#!/usr/bin/env python3
import argparse, sys, os, re, logging, traceback
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def protein_dict(fasta_prot):
    """
    Filters proteins by length and replaces ambiguous amino acids.
    
    Args:
        fasta_prot (str): Path to input fasta file.
    Return:
        dict: Filtered proteins with updated sequences.
    """
    prots = dict()
    # Use BioIPython SeqIO to parse the fasta
    with open(fasta_prot, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # 1. Length filter: skip proteins shorter than 30 aa
            if len(record.seq) < 30:
                continue
            
            # 2. Sequence cleaning: 
            # Replace ambiguous AA (U, O, B, Z) with '@'
            clean_seq = re.sub(r"[UOBZ]", "@", str(record.seq))
            
            # Store in dictionary
            prots[record.id] = clean_seq
            
    return prots

def main():
    try:
        parser = argparse.ArgumentParser(description="NLS Filter and Pre-processing")
        parser.add_argument("--input", required=True, help="Input protein fasta file")
        parser.add_argument("--output", required=True, help="Output filtered fasta file")
        args = parser.parse_args()

        prot_dict = protein_dict(args.input)
        
        # Write output 
        filtered_records = []
        for p_id, p_seq in prot_dict.items():
            record = SeqRecord(Seq(p_seq), id=p_id, description="")
            filtered_records.append(record)
        
        with open(args.output, "w") as output_handle:
            SeqIO.write(filtered_records, output_handle, "fasta")
            
        logging.info(f"Filtering complete. {len(filtered_records)} proteins saved to {args.output}")

    except Exception as ex:
        logging.error(f"Error during NLS filtering: {str(ex)}")
        traceback.print_exc()
        sys.exit(os.EX_SOFTWARE)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
