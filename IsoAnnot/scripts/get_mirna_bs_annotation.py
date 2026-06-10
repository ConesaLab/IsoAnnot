#!/usr/bin/env python3
"""
Author: Raúl López

Script to map genomic miRNA binding sites (from miRWalk) to transcriptomic coordinates. 
It performs a sequence-level validation by comparing the genomic sequence provided 
by miRWalk with the actual sequence extracted from the transcript FASTA file, 
ensuring the mapping is accurate across exons and introns.
"""

import argparse, sys, csv, logging, traceback
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict
from IsoAnnot import read_chr_ref_acc
from bx.intervals.intersection import IntervalTree

def load_fasta_seqs(fasta_file):
    """
    Loads FASTA sequences and creates two indices to handle naming discrepancies 
    between annotation versions (e.g., presence or absence of version suffixes).

    Args:
        fasta_file (str): Path to the isoform FASTA file.

    Returns:
        tuple: (seqs_exact, seqs_no_ver) dictionaries where keys are IDs and values are sequences.
    """
    logging.info(f"Loading sequence of isoforms: {fasta_file}")
    seqs_exact = {}
    seqs_no_ver = {}

    with open(fasta_file) as handle:
        for record in SimpleFastaParser(handle):
            # Extract full ID (stripping description after first whitespace)
            full_id = record[0].split()[0]
            sequence = record[1]

            # Store exact match (e.g., ENST000001.1)
            seqs_exact[full_id] = sequence

            # Store version-less match (e.g., ENST000001) to increase compatibility
            if '.' in full_id:
                clean_id = full_id.split('.')[0]
                seqs_no_ver[clean_id] = sequence

    return seqs_exact, seqs_no_ver

def load_genepred(file_path, chr_ref=None):
    """
    Loads transcript structures from a GenePred file into an IntervalTree for 
    fast genomic overlap lookups.

    Args:
        file_path (str): Path to the GenePred file.
        chr_ref (dict, optional): Dictionary for chromosome name mapping (NCBI to Ensembl).

    Returns:
        tuple: (isoforms_by_chr, exons_by_chr) dictionaries for structural analysis.
    """
    logging.info(f"Loading structure (GenePred): {file_path}")
    data = defaultdict(dict)
    exons_by_chr = defaultdict(lambda: IntervalTree())

    if chr_ref is None:
        chr_ref = {}

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            cols = line.strip().split("\t")
            t_id = cols[0]
            chrom = cols[1]
            # Map chromosome name if conversion table is provided
            clean_chrom = chr_ref.get(chrom, chrom)
            strand = cols[2]
            starts = [int(x) for x in cols[8].strip(",").split(",")]
            ends = [int(x) for x in cols[9].strip(",").split(",")]
            exons = list(zip(starts, ends))
            
            # Order exons based on strand for transcriptomic coordinate calculation
            if strand == "+":
                biographical_exons = exons
            else:
                biographical_exons = exons[::-1]

            data[clean_chrom][t_id] = {'exons': biographical_exons, 'strand': strand}
            # Insert every exon into the Interval Tree for the specific chromosome
            for s, e in exons:
                exons_by_chr[clean_chrom].insert(s, e, t_id)
    return data, exons_by_chr

def get_transcript_length_at_genomic_pos(ordered_exons, genomic_pos, strand):
    """
    Translates a 1-based genomic position into a 1-based transcriptomic coordinate, 
    accounting for introns and strand orientation.

    Args:
        ordered_exons (list): List of (start, end) tuples ordered by transcription flow.
        genomic_pos (int): The genomic coordinate to translate.
        strand (str): Strand of the transcript ('+' or '-').

    Returns:
        int: The relative position within the transcript, or -1 if the position is intronic.
    """
    transcript_pos = 0

    for exon_start, exon_end in ordered_exons:
        exon_len = exon_end - exon_start
        # Check if genomic position falls within the current exon
        if exon_start <= (genomic_pos - 1) < exon_end:
            if strand == "+":
                offset = genomic_pos - exon_start
            else:
                offset = exon_end - genomic_pos + 1
            return transcript_pos + offset

        # Add current exon length to the running total
        transcript_pos += exon_len
    return -1

def main():
    """
    Main execution flow:
    1. Load genomic structures and sequences.
    2. Parse miRWalk genomic predictions.
    3. Project genomic coordinates to transcripts.
    4. Cross-validate sequences against FASTA records.
    5. Write validated binding sites to a GTF-like output.
    """
    parser = argparse.ArgumentParser(description="Map genomic miRNA BS to transcriptomic coordinates")
    parser.add_argument("--genepred", required=True, help="GenePred file of isoforms")
    parser.add_argument("--isoform_fasta", required=True, help="FASTA file of isoforms")
    parser.add_argument("--mirwalk_genomic", required=True, help="miRWalk file with genomic coordinates")
    parser.add_argument("--output", required=True, help="Final output TSV/GTF")
    parser.add_argument("--db", nargs=None, required=True, help="Database type (ensembl/refseq)")
    parser.add_argument("--chr_ref", required=False, help="Refseq-Ensembl chromosome mapping table")
    args = parser.parse_args()

    try:
        # Step 1: Handle chromosome accessions mapping for RefSeq
        refseqChrom = {}
        if args.chr_ref and args.db.lower() == "refseq":
            refseqChrom = read_chr_ref_acc(args.chr_ref)
            logging.info("Using Refseq chromosome accessions mapping.")
            isoforms_by_chr, exons_by_chr = load_genepred(args.genepred, refseqChrom)
        else:
            logging.info("Using original chromosome accessions.")
            isoforms_by_chr, exons_by_chr = load_genepred(args.genepred)

        # Step 2: Load sequences with double indexing (Exact vs No-Version)
        isoform_seqs_exact, isoform_seqs_no_ver = load_fasta_seqs(args.isoform_fasta)

        logging.info("Starting sequence validation process...")

        with open(args.mirwalk_genomic, "r") as f_in, open(args.output, "w") as f_out:
            reader = csv.DictReader(f_in, delimiter="\t")
            writer = csv.writer(f_out, delimiter="\t")
            count_valid = 0
            count_mismatch = 0
            
            for row in reader:
                overlaps = set()
                miRNA = row["miRNA"]
                chrom = row["Chrom"]
                g_start = int(row["G_Start"])
                g_end = int(row["G_End"])
                strand = row["Strand"]
                source = row["Source"]
                expected_sequence = row["Sequence"]

                # Chromosome sanity check
                if chrom not in isoforms_by_chr or chrom not in exons_by_chr: 
                    continue

                # Find all isoforms overlapping the genomic range of the binding site
                for interval in exons_by_chr[chrom].find(g_start, g_end):
                    overlaps.add(interval)
                
                for iso_id in overlaps:
                    iso_data = isoforms_by_chr[chrom][iso_id]
                    
                    # Only map if strands are consistent
                    if iso_data['strand'] == strand:
                        ex_strand = iso_data['strand']
                        exons = iso_data['exons']
                        
                        # Project genomic start/end to transcript coordinates
                        t_start = get_transcript_length_at_genomic_pos(exons, g_start, ex_strand)
                        t_end = get_transcript_length_at_genomic_pos(exons, g_end, ex_strand)

                        if t_start != -1 and t_end != -1:
                            final_t_start = min(t_start, t_end)
                            final_t_end = max(t_start, t_end)
                            
                            actual_sequence = ""
                            full_seq = ""

                            # Handle naming mismatches with 3 fallback attempts
                            # Attempt 1: Exact match (ENST001.1)
                            if iso_id in isoform_seqs_exact:
                                full_seq = isoform_seqs_exact[iso_id]
                                actual_sequence = full_seq[final_t_start-1 : final_t_end]

                            # Attempt 2: GenePred has version, FASTA ID does not
                            elif iso_id.split('.')[0] in isoform_seqs_exact:
                                full_seq = isoform_seqs_exact[iso_id.split('.')[0]]
                                actual_sequence = full_seq[final_t_start-1 : final_t_end]

                            # Attempt 3: GenePred has no version, FASTA ID has
                            elif iso_id in isoform_seqs_no_ver:
                                full_seq = isoform_seqs_no_ver[iso_id]
                                actual_sequence = full_seq[final_t_start-1 : final_t_end]

                            # Sequence validation (Case-insensitive)
                            if full_seq:
                                if actual_sequence.upper() == expected_sequence.upper():
                                    count_valid += 1
                                    gtf_attributes = f"ID={miRNA}"
                                    writer.writerow([
                                        iso_id, source, "miRNA_binding_site", 
                                        final_t_start, final_t_end,
                                        ".", ex_strand , ".", gtf_attributes
                                    ])
                                else:
                                    count_mismatch += 1
                                    # Useful for debugging coordinate shift issues
                                    logging.warning(f"Mismatch in {iso_id}: Expected {expected_sequence} but got {actual_sequence}")

        logging.info(f"Validation completed. Valid: {count_valid}, Discarded: {count_mismatch}")

    except Exception as ex:
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    # Clean logging format for Snakemake logs
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
