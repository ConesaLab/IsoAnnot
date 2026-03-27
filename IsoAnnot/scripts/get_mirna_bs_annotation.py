#!/usr/bin/env python3
import argparse, sys, csv, logging, traceback
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict
from IsoAnnot import read_chr_ref_acc
from bx.intervals.intersection import IntervalTree

def load_fasta_seqs(fasta_file):
    """
    Loads FASTA sequences and creates two indices for maximum flexibility:
    1. seqs_exact: Full ID -> Sequence (e.g., ENST000.1 -> Seq)
    2. seqs_no_ver: ID without version -> Sequence (e.g., ENST000 -> Seq)
    """
    logging.info(f"Loading sequence of isoforms: {fasta_file}")
    seqs_exact = {}
    seqs_no_ver = {}
    
    with open(fasta_file) as handle:
        for record in SimpleFastaParser(handle):
            # Full ID as it appears in the FASTA (removing description after whitespace)
            full_id = record[0].split()[0]
            sequence = record[1]
            
            # 1. Store exact match
            seqs_exact[full_id] = sequence
            
            # 2. Store match without version (for Attempt 3 in main loop)
            # Only if it has a dot, we store the "clean" version pointing to the sequence.
            # NOTE: For PacBio (PB.1.1), this maps 'PB' -> Seq. 
            # This is safe because PacBio will always match in Step 1 (Exact).
            if '.' in full_id:
                clean_id = full_id.split('.')[0]
                seqs_no_ver[clean_id] = sequence
                
    return seqs_exact, seqs_no_ver

def load_genepred(file_path, chr_ref=None):
    """
    Loads GenePred structure. 
    Accepts an optional 'chr_ref' dictionary for chromosome mapping.
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
            clean_chrom = chr_ref.get(chrom, chrom)
            strand = cols[2]
            starts = [int(x) for x in cols[8].strip(",").split(",")]
            ends = [int(x) for x in cols[9].strip(",").split(",")]
            exons = list(zip(starts, ends))
            if strand == "+":
                biographical_exons = exons
            else:
                biographical_exons = exons[::-1]

            data[clean_chrom][t_id] = {'exons': biographical_exons, 'strand': strand}
            for s, e in exons:
                exons_by_chr[clean_chrom].insert(s, e, t_id)
    return data, exons_by_chr


def get_transcript_length_at_genomic_pos(ordered_exons, genomic_pos, strand):
    """
    Calculates transcriptomic coordinates (1-based).
    """
    transcript_pos = 0
    
    for exon_start, exon_end in ordered_exons:
        exon_len = exon_end - exon_start
        if exon_start <= (genomic_pos - 1) < exon_end:
            if strand == "+":
                offset = genomic_pos - exon_start
            else:
                offset = exon_end - genomic_pos + 1
            return transcript_pos + offset
            
        transcript_pos += exon_len
    return -1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genepred", required=True, help="GenePred file of your isoforms")
    parser.add_argument("--isoform_fasta", required=True, help="FASTA file of your isoforms")
    parser.add_argument("--mirwalk_genomic", required=True, help="miRWalk file filtered by species at genomic coordinates")
    parser.add_argument("--output", required=True, help="Final output")
    parser.add_argument("--db", nargs=None, required=True)
    parser.add_argument("--chr_ref", required=False, help=" Refseq chromosome accessions and their corresponding Ensembl identifyer")
    args = parser.parse_args()

    try:
        refseqChrom = {}
        if args.chr_ref and args.db.lower() == "refseq":
            refseqChrom = read_chr_ref_acc(args.chr_ref)
            logging.info("Using Refseq chromosome accessions mapping.")
            isoforms_by_chr, exons_by_chr = load_genepred(args.genepred, refseqChrom)
        else:
            logging.info("Using original chromosome accessions.")
            isoforms_by_chr, exons_by_chr = load_genepred(args.genepred)
        
        isoform_seqs_exact, isoform_seqs_no_ver = load_fasta_seqs(args.isoform_fasta)
        
        logging.info("Starting validation of sequence...")



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
                
                if chrom not in isoforms_by_chr or chrom not in exons_by_chr: continue
                
                if len(exons_by_chr[chrom].find(g_start, g_end)) == 0:
                    continue
                for interval in exons_by_chr[chrom].find(g_start, g_end):
                    overlaps.add(interval)
                for iso_id in overlaps:
                    iso_data = isoforms_by_chr[chrom][iso_id]
                    if iso_data['strand'] == strand:
                        ex_strand = iso_data['strand']
                        exons = iso_data['exons']
                        t_start = get_transcript_length_at_genomic_pos(exons, g_start, ex_strand)
                        t_end = get_transcript_length_at_genomic_pos(exons, g_end, ex_strand)
                        
                        if t_start != -1 and t_end != -1:
                            final_t_start = min(t_start, t_end)
                            final_t_end = max(t_start, t_end)
                
                            # Attempt 1: Exact Match 
                            if iso_id in isoform_seqs_exact:
                                full_seq = isoform_seqs_exact[iso_id]
                                actual_sequence = full_seq[final_t_start-1 : final_t_end]
                            
                            # Attempt 2: GP has version, FASTA does not
                            elif iso_id.split('.')[0] in isoform_seqs_exact:
                                full_seq = isoform_seqs_exact[iso_id.split('.')[0]]
                                actual_sequence = full_seq[final_t_start-1 : final_t_end]
                                
                            # Attempt 3: GP has no version, FASTA has
                            # Use the auxiliary 'isoform_seqs_no_ver' dictionary
                            elif iso_id in isoform_seqs_no_ver:
                                full_seq = isoform_seqs_no_ver[iso_id]
                                actual_sequence = full_seq[final_t_start-1 : final_t_end]

                            if full_seq:
                                if actual_sequence.upper() == expected_sequence.upper():
                                    count_valid += 1
                                    gtf_attributes = f"ID={miRNA}"
                                    writer.writerow([
                                        iso_id, source, "miRNA_binding_site", final_t_start, final_t_end,
                                        ".", ex_strand , ".", gtf_attributes
                                    ])
                                else:
                                    count_mismatch += 1
                                    logging.warning(f"\n--- DEBUG MISMATCH ---")
                                    logging.warning(f"Isoforma: {iso_id} | Strand: {ex_strand}")
                                    logging.warning(f"Genomic: {g_start}-{g_end} | Calc Transcript: {final_t_start}-{final_t_end}")
                                    logging.warning(f"Expected (miRWalk): {expected_sequence}")
                                    logging.warning(f"Actual   (FASTA)  : {actual_sequence}")

        logging.info(f"Validation completed.")
        logging.info(f"miRNA binding sites validated: {count_valid}")
        logging.info(f"miRNA binding sites discarded: {count_mismatch}")

    except Exception as ex:
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
