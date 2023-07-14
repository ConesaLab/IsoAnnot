#!/usr/bin/env python3
'''
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez
Modified by: Alessandra Martinez

Script that gets the NMD status and the reference
proteins associated to each sequenced transcript.
'''

from collections import defaultdict
import argparse, sys, os, csv, logging, re, traceback, json
import pandas as pd
from IsoAnnot import strand_table, get_summary_parts, read_chr_ref_acc
from Bio.SeqFeature import FeatureLocation, SeqFeature
from uniprotPhosphosite_genomicCoordinates import get_fasta_sequences
from layer_exons import get_gtf_exons


def _get_transcript_exons_location(gtf_data, transcript):
    '''
    Constructs the location of the exons

    Args:
        gtf_data ():
        transcript_id (str)

    Returns:
        SeqFeature for exon location

    '''
    isoform_exons = []
    for line in gtf_data:
        transcript_id = line.split("transcript_id")[1].split(";")[0].replace('"', '').strip()  
        if transcript_id == transcript:
            chr = line.split("\t")[0]
            transcript_exon_start=int(line.split("\t")[3])
            transcript_exon_end=int(line.split("\t")[4])
            transcript_exon_strand = line.split("\t")[6]
            isoform_exons.append({"start" :transcript_exon_start, "end": transcript_exon_end, "strand": transcript_exon_strand})

    isoforms_exon_locs = [
        FeatureLocation(exon_block['start'], exon_block['end'], strand_table.get(exon_block['strand'], None)) for exon_block
        in isoform_exons
    ]

    if len(isoforms_exon_locs) > 1:
        final_location = sum(isoforms_exon_locs)
    else:
        final_location = isoforms_exon_locs.pop()

    return SeqFeature(final_location, id=chr)


def _calculate_nmd(prot_start, prot_end, exon_locations):
    '''
    Calculates nonsense mediated decay with the 55 nucleotide rule:
    A termination codon that falls more than 55 nt upstream of an exon-exon junction 
    will be potentially detected and degraded via the NMD machinery.
    '''
    exons_starts, exons_ends = [], []

    for part in exon_locations.parts:
        exons_starts.append(part.start)
        exons_ends.append(part.end)

    exons_starts.sort()
    exons_ends.sort()

    prot_ord = sorted([prot_start, prot_end])

    nmd_status = "nonNMD"

    if exon_locations.strand == 1 and prot_ord[1] < exons_starts[-1] and (exons_ends[-2] - prot_ord[1] + 1) > 55:
        nmd_status = "NMD"
    elif exon_locations.strand == -1 and prot_ord[0] > exons_ends[0] and (prot_ord[0] - exons_starts[1] + 1) > 55:
        nmd_status = "NMD"

    return nmd_status


def _calculate_cds_positions(cds_start, cds_end, location):
    '''
    Gets CDS positions and NMD status for a transcript
    '''
    exon_sorted_positions = sorted(location, reverse=location.strand < 0)
    cds_positions = sorted(exon_sorted_positions[cds_start - 1:cds_end])  # Ver si esta bien el restar 1

    output_dict = get_summary_parts(cds_positions)
    output_dict["nmd"] = _calculate_nmd(cds_positions[0], cds_positions[-1], location)

    return output_dict


def _retrieve_pacbio_info(gtf_data, fasta_data):  
    '''
    Gets header information from predicted ORFs (isoform ID and CDS)
    '''
    capture_regex = re.compile(r"(.+)\s.*?\|.*?\|.*?\|.*?\|([0-9]+)")
    isoforms_cds = {}

    for fasta_id, fasta_seq in fasta_data.items():
        match_id = re.match(capture_regex, fasta_id)

        if match_id:
            isoform_id, cds_start = match_id.groups()
            cds_start = int(cds_start)
            cds_end = cds_start + len(fasta_seq)*3 -1 
            exons_final_location = _get_transcript_exons_location(gtf_data, isoform_id)
            isoforms_cds[isoform_id] = _calculate_cds_positions(cds_start, cds_end, exons_final_location.location)
            isoforms_cds[isoform_id]['chr'] = exons_final_location.id

    return isoforms_cds


def _retrieve_pacbio_id(fasta_data):
    '''
    Cleans a fasta dictionary by keeping only the isoform
    ID from the headers
    '''
    capture_regex = re.compile(r"(.+)\s.*")

    def _id(old_id):
        return re.match(capture_regex, old_id).group(1)
    
    cleaned_dict = {_id(k): v for k, v in fasta_data.items()}

    return cleaned_dict


def read_genomic_cds_from_gtf_group_chr(gtf_file, chr_ref={}):
    '''
    Gets CDS genomic positions 
    '''
    aux_dict = {"Chromosome":[], "protein_id": [], "Start":[], "End": []}
    parsed_dict = defaultdict(dict)

    with open(gtf_file, "r") as gtf_content:
        for line in gtf_content:
            if line[0] != "#":
                if line.split("\t")[2]=="CDS":
                    try:
                        protein_id = line.split("protein_id")[1].split(";")[0].replace('"', '').strip()
                    except IndexError:
                        continue
                    chr = line.split("\t")[0]
                    chr_name = chr_ref.get(chr, chr)
                    start = line.split("\t")[3]
                    end = line.split("\t")[4]
                    aux_dict["Chromosome"].append(chr_name)
                    aux_dict["protein_id"].append(protein_id)
                    aux_dict["Start"].append(start)
                    aux_dict["End"].append(end)
    
    gtf_cds = pd.DataFrame.from_dict(aux_dict)

    for chr, values in gtf_cds.groupby("Chromosome"):
        for prot_id, prot_values in values.groupby("protein_id"):
            parsed_dict[chr][prot_id] = {
                "start": sorted(prot_values['Start'].tolist()),
                "end": sorted(prot_values['End'].tolist())
            }

    return parsed_dict


def main():

    # Inputs
    parser = argparse.ArgumentParser(description="PacbioToReference")
    parser.add_argument("--ensembl_gtf", required=True)
    parser.add_argument("--refseq_gtf", required=False, default=None)
    parser.add_argument("--chr_ref", required=False, default=None)

    parser.add_argument("--classification_file", required=False)
    parser.add_argument("--corrected_gtf", required=True)
    parser.add_argument("--corrected_fasta_proteins", required=True)

    parser.add_argument("--database", required=True)
    parser.add_argument("--species_db", required=True)

    # Outputs
    parser.add_argument("--output_assoc", required=True, type=str)
    parser.add_argument("--output_nmd", required=True, type=str)

    args = parser.parse_args()

    try: 
        
        # Read reference proteins information (Uniprot, Ensembl, RefSeq)
        with open(args.species_db, "rb") as fp:
            all_proteins = json.load(fp)

        # Get chromosome Ensembl-RefSeq IDs equivalence if we
        # are using GTF data
        if args.refseq_gtf:
            chrom_ref = read_chr_ref_acc(args.chr_ref)

        # Get Ensembl and RefSeq CDS info
        cds_ensembl = read_genomic_cds_from_gtf_group_chr(args.ensembl_gtf)
        if args.refseq_gtf:
            cds_refseq = read_genomic_cds_from_gtf_group_chr(args.refseq_gtf,
                                      chrom_ref)
            
        # Get isoforms exon locations and predicted proteins of sequenced transcriptome
        # using SQANTI outputs
        corrected_exons = get_gtf_exons(args.corrected_gtf)
        corrected_fasta_proteins = get_fasta_sequences(fasta_files=[args.corrected_fasta_proteins],
                                               matching_regex="(.*)",
                                               keep_version="True")

        # Get cds positions of the transcripts
        cds_nmd_info = _retrieve_pacbio_info(corrected_exons, corrected_fasta_proteins)

        # In the previous step we needed the original fasta headers to extract useful info (ie. cds positions
        # for Pacbio). Now we clean the data by keeping only the transcript id.
        corrected_fasta_proteins = _retrieve_pacbio_id(corrected_fasta_proteins)

        # 2 output files: associated proteins and NMD
        with open(args.output_assoc, "w") as output_file_assoc, open(args.output_nmd, "w") as output_file_nmd:
            output_assoc_tsv = csv.writer(output_file_assoc, delimiter="\t")
            output_nmd_tsv = csv.writer(output_file_nmd, delimiter="\t")

            # For coding isoforms
            for isoform_id, isoform_info in cds_nmd_info.items():
                if isoform_id in corrected_fasta_proteins:

                    isoform_seq = corrected_fasta_proteins.get(isoform_id)
                    associated_prots = []

                    # Get nmd info
                    nmd_status = isoform_info.get("nmd")
                    output_nmd_tsv.writerow([isoform_id, nmd_status])

                    # Keep reference proteins that have the same sequence as the predicted isoform protein 
                    for source_tag, protein_ids in all_proteins.get(isoform_seq, {}).items():
                        for protein in protein_ids:
                            associated_prots.append(protein)
                            output_assoc_tsv.writerow([isoform_id, protein, source_tag, "sequence"])

                    # # Compare CDSs between Pacbio and reference transcripts and keep proteins that match CDS coordinates
                    # transcript_starts = isoform_info.get('start')
                    # transcript_ends = isoform_info.get('end')
                    # transcript_chr = isoform_info.get('chr')

                    # chr_ensembl = cds_ensembl.get(transcript_chr, {})
                    # for cds_ref, cds_props in chr_ensembl.items():
                    #     if cds_props.get("start") == transcript_starts and cds_props.get("end") == transcript_ends:
                    #         output_assoc_tsv.writerow([isoform_id, cds_ref, "ensembl", "coord"])

                    # if cds_refseq:
                    #     chr_refseq = cds_refseq.get(transcript_chr, {})
                    #     for cds_ref, cds_props in chr_refseq.items():
                    #         if cds_props.get("start") == transcript_starts and cds_props.get("end") == transcript_ends:
                    #             output_assoc_tsv.writerow([isoform_id, cds_ref, "refseq", "coord"])

    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        traceback.print_exc()
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)

if __name__ == "__main__":
    main()
