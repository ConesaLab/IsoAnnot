#!/usr/bin/env python3
"""
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez
Modified: Alessandra Martínez

Script to get exons structural information and store it in GTF format.
"""

import argparse, sys, os, logging, csv
from IsoAnnot import read_chr_ref_acc, argparse_nullable


def get_gtf_exons(gtf_file:str):
    """
    Filters lines of GTF containing exon information.

    Args: 
        gtf_file (str): GTF file

    Returns:
        gtf_exons (list): of .
    """
    gtf_exons = []
    with open(gtf_file, "r") as gtf_content:
        for line in gtf_content:
            if line[0] != "#":
                if line.split("\t")[2]=="exon":
                    try:
                        transcript_id = line.split("transcript_id")[1].split(";")[0].replace('"', '').strip()
                        gtf_exons.append(line)
                    except IndexError:
                        continue
    return gtf_exons
        

def main():
    try:
        parser = argparse.ArgumentParser(description='IsoAnnot Exons Layer')
        parser.add_argument('--gtf_file', nargs=None, required=True)
        parser.add_argument('--output', nargs=None, required=True)
        parser.add_argument('--chr_ref', required=False, nargs='?', const=None, type=argparse_nullable)

        args = parser.parse_args()
        
        # Get info from file
        exons_info = get_gtf_exons(args.gtf_file)
        # exons_info = get_gtf_exon_locations(args.gtf_file)

        # Get chromosome Ensembl-Refseq IDs
        chrom_ref = {}
        if args.chr_ref:
            chrom_ref=read_chr_ref_acc(args.chr_ref)

        with open(args.output, "w") as output_file:
            logging.info(f"Starting exons conversion to GTF {args.output}")
            output_file_tsv = csv.writer(output_file, delimiter="\t")

            for line in exons_info: 
                chr = line.split("\t")[0]
                chr_name = chrom_ref.get(chr, chr)
                transcript_id = line.split("transcript_id")[1].split(";")[0].replace('"', '').strip()
                transcript_exon_start=int(line.split("\t")[3])
                transcript_exon_end=int(line.split("\t")[4])
                transcript_exon_strand = line.split("\t")[6]
                output_file_tsv.writerow([
                        transcript_id,  # seqname
                        "tappAS",  # source
                        "exon",  # feature
                        transcript_exon_start,  # start
                        transcript_exon_end,  # end
                        ".",  # score
                        transcript_exon_strand,  # strand
                        ".",  # frame
                        f"Chr=chr{chr_name}"  # attribute
                    ])

    except Exception as ex:
        logging.error(str(ex))
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
