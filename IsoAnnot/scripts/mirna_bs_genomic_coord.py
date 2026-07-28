#!/usr/bin/env python3
"""
Author: Raúl López

Script to map transcriptomic miRNA binding sites (from miRWalk) to genomic coordinates.
It uses a reference GTF to understand exon structures and a FASTA file to validate 
that the transcriptomic sequence at the predicted binding site matches the miRWalk record.
"""

import argparse, sys, csv, logging, traceback, re
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict
from IsoAnnot import read_chr_ref_acc

# Standard logging configuration for pipeline traceability
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def get_exon_from_gtf(gtf_file, chrom_ref={}):
    """
    Parses a GTF file to extract exon coordinates grouped by transcript.

    Args:
        gtf_file (str): Path to the reference GTF file.
        chrom_ref (dict): Dictionary for chromosome name conversion (e.g., NCBI to Ensembl).

    Returns:
        dict: A dictionary where keys are transcript IDs and values are lists of 
              sorted exon coordinate dictionaries.
    """
    logging.info(f"Reading GTF: {gtf_file}...")
    transcript_coord = defaultdict(list)
    ids_loaded = 0

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if parts[2] != "exon": continue

            attributes = parts[8]
            try:
                # Extract transcript_id using regex
                t_match = re.search(r'transcript_id "([^"]+)"', attributes)
                if not t_match: continue
                t_id_raw = t_match.group(1)
            except AttributeError:
                continue

            # Standardize ID (remove version) and handle chromosome mapping
            t_id_raw = t_id_raw.split(".")[0]
            raw_chrom = parts[0]
            clean_chrom = chrom_ref.get(raw_chrom, raw_chrom)

            seq_record = {
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "chr": clean_chrom
            }
            transcript_coord[t_id_raw].append(seq_record)
            ids_loaded += 1

    # Ensure exons are sorted by genomic start position for consistent mapping
    for t_id in transcript_coord:
        transcript_coord[t_id].sort(key=lambda x: x["start"])

    logging.info(f"GTF loaded. {len(transcript_coord)} unique transcripts indexed.")
    return transcript_coord

def get_fasta_seq(refseq_fasta):
    """
    Loads sequences from a FASTA file into a dictionary.

    Args:
        refseq_fasta (str): Path to the FASTA file.

    Returns:
        dict: Keys are cleaned sequence IDs, values are the sequences strings.
    """
    logging.info(f"Reading FASTA: {refseq_fasta}...")
    output_results = dict()
    capture_regex = re.compile(r">?([\w\.]+)")
    with open(refseq_fasta) as handle:
        for record in SimpleFastaParser(handle):
            match = capture_regex.match(record[0])
            if match:
                # Clean version from ID for matching
                id_name = match.group(1).split(".")[0]
                output_results[id_name] = record[1]

    return output_results

def map_transcript_to_genomic_blocks(transcript_exons, t_start, t_end):
    """
    Translates transcriptomic coordinates (linear) back to genomic blocks (discontinuous).
    This accounts for splicing events and strand orientation.

    Args:
        transcript_exons (list): Sorted list of exon dictionaries for a transcript.
        t_start (int): 1-based start position relative to the transcript.
        t_end (int): 1-based end position relative to the transcript.

    Returns:
        list: Sorted list of (start, end) tuples representing genomic segments of the binding site.
    """
    transcript_exons.sort(key=lambda x: x['start'])
    strand = transcript_exons[0]['strand']

    # Biological order: exons must be processed 5' to 3'
    if strand == "-":
        biologic_exons = transcript_exons[::-1]
    else:
        biologic_exons = transcript_exons

    genomic_blocks = []
    current_transcript_pos = 0

    for exon in biologic_exons:
        exon_len = (exon["end"] - exon["start"]) + 1
        exon_t_start = current_transcript_pos + 1
        exon_t_end = current_transcript_pos + exon_len

        # Determine if the binding site range overlaps with the current exon
        overlap_start = max(t_start, exon_t_start)
        overlap_end = min(t_end, exon_t_end)

        if overlap_start <= overlap_end:
            # Calculate distance from the exon's transcriptomic start
            offset_start = overlap_start - exon_t_start
            offset_end = overlap_end - exon_t_start

            # Translate offsets to genomic positions based on strand
            if strand == "+":
                g_block_start = exon["start"] + offset_start
                g_block_end = exon["start"] + offset_end
            else:
                g_block_start = exon["end"] - offset_end
                g_block_end = exon["end"] - offset_start

            genomic_blocks.append((g_block_start, g_block_end))

        current_transcript_pos += exon_len

    # Final genomic blocks must be sorted by start coordinate for output
    genomic_blocks.sort(key=lambda x: x[0])
    return genomic_blocks

def main():
    """
    Main execution logic:
    1. Parse miRWalk interaction file.
    2. Retrieve structural and sequence info for mRNA targets.
    3. Validate binding site sequences.
    4. Project transcript coordinates to genomic ranges.
    5. Output a comprehensive mapping TSV.
    """
    parser = argparse.ArgumentParser(description="Map transcriptomic miRWalk interactions to genomic coordinates")
    parser.add_argument("--mirwalk_file", required=True, help="Filtered miRWalk output")
    parser.add_argument("--refseq_fasta", required=True, help="Reference FASTA for targets")
    parser.add_argument("--refseq_gtf", required=True, help="Reference GTF for targets")
    parser.add_argument("--mirna_output", required=True, help="Output path for genomic mapping")
    parser.add_argument("--chr_ref", required=False, help="Chromosome accession mapping file")
    args = parser.parse_args()

    try:
        # Load chromosome mapping if provided (e.g., NCBI accessions to Ensembl names)
        refseqChrom = read_chr_ref_acc(args.chr_ref, leading_db="ensembl") if args.chr_ref else {}

        coord_dict = get_exon_from_gtf(args.refseq_gtf, chrom_ref=refseqChrom)
        
        if hasattr(refseqChrom, 'total_queries') and refseqChrom.total_queries > 0:
            failure_rate = len(refseqChrom.unmapped_log) / refseqChrom.total_queries
            if failure_rate > 0.50:
                raise RuntimeError(
                    f"Critical Error in mirna_bs_genomic_coord.py: {len(refseqChrom.unmapped_log)}/{refseqChrom.total_queries} "
                    f"chromosome lookups failed using mapping file '{args.chr_ref}'. Halting execution to prevent silent data loss."
                )

        fasta_dict = get_fasta_seq(args.refseq_fasta)

        logging.info("Starting coordinate mapping and sequence validation...")

        count_mapped = 0
        count_total = 0
        count_no_gtf = 0
        count_no_fasta = 0
        count_seq_mismatch = 0
        count_pos_mismatch = 0

        with open(args.mirna_output, "w", newline='') as f_out, open(args.mirwalk_file, "r") as f_in:
            writer = csv.writer(f_out, delimiter="\t")

            header = ["FeatureType", "miRNA", "mRNA", "Chrom", "Strand", "G_Start", "G_End",
                      "Exon_Starts", "Exon_Ends", "Len_NT", "Sequence", "Source"]
            writer.writerow(header)

            for line in f_in:
                if line.startswith("miRNA"): continue
                cols = line.strip().split("\t")
                if len(cols) < 4: continue

                count_total += 1
                mirna = cols[0]
                t_id = cols[1]

                # Basic validation: IDs must exist in reference files
                if t_id not in coord_dict:
                    count_no_gtf += 1
                    continue

                # Extract the expected duplex sequence from miRWalk field
                duplex_field = cols[4]
                if "#" in duplex_field:
                    expected_sequence = duplex_field.split("#")[1]
                else:
                    continue
                
                motif = "."
                sequence_match = False
                
                if t_id in fasta_dict:
                    try:
                        # Extract transcriptomic coordinates
                        bs_raw_col = cols[3]
                        bs_raw = bs_raw_col.split(",")
                        bs_start, bs_end = int(bs_raw[0]), int(bs_raw[1])
                        
                        # Extract sequence from FASTA to cross-validate
                        motif = fasta_dict[t_id][bs_start-1:bs_end]
                        
                        # NOTE: miRWalk sequences often skip the first nucleotide 
                        # of the mapping, hence the [1:] comparison.
                        if motif.upper()[1:] == expected_sequence.upper():
                            sequence_match = True
                    except Exception as e:
                        logging.warning(f"Error extracting sequence: {e}")
                        pass
                else:
                    count_no_fasta += 1
                    continue

                if not sequence_match:
                    count_seq_mismatch += 1
                    continue

                try:
                    # Sanity check: coordinates must not exceed total transcript length
                    total_transcript_len = sum([(exon["end"] - exon["start"] + 1) for exon in coord_dict[t_id]])

                    if bs_end > total_transcript_len:
                        count_pos_mismatch += 1
                        continue

                    # Translate linear binding site to genomic blocks
                    g_blocks = map_transcript_to_genomic_blocks(coord_dict[t_id], bs_start, bs_end)

                    if g_blocks:
                        exon_info = coord_dict[t_id][0]
                        chrom = exon_info["chr"]
                        strand = exon_info["strand"]

                        global_start = g_blocks[0][0]
                        global_end = g_blocks[-1][1]
                        starts_str = ",".join([str(b[0]) for b in g_blocks])
                        ends_str = ",".join([str(b[1]) for b in g_blocks])
                        length = bs_end - bs_start + 1

                        writer.writerow([
                            "miRNA_Binding_Site", mirna, t_id, chrom, strand,
                            global_start, global_end, starts_str, ends_str,
                            length, motif, "miRWalk"
                        ])
                        count_mapped += 1
                except Exception as e:
                    logging.warning(f"Failed to map binding site for transcript {t_id}: {e}")
                    continue

            # Final summary statistics
            logging.info("-" * 30)
            logging.info(f"FINAL RESULT SUMMARY:")
            logging.info(f"Binding sites read: {count_total}")
            logging.info(f"Binding sites correctly mapped: {count_mapped}")
            logging.info(f"Failed: Missing GTF ID: {count_no_gtf}")
            logging.info(f"Failed: Missing FASTA ID: {count_no_fasta}")
            logging.info(f"Failed: Sequence mismatch: {count_seq_mismatch}")
            logging.info(f"Failed: Index out of range: {count_pos_mismatch}")

    except Exception as ex:
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
