#!/usr/bin/env python3
# Author: Lorena de la Fuente Lorente
# Refactored: Carlos Martínez
# Modified by Fran Pardo-Palacios and Alessandra Martinez

import argparse, csv, logging, sys, os, traceback, json
from IsoAnnot import merge_fasta_dicts
from uniprotPhosphosite_genomicCoordinates import get_fasta_sequences


def main():
    try: 
        parser = argparse.ArgumentParser(description="Parse Databases")
        parser.add_argument("--uniprot_fasta", nargs="+", required=True)
        parser.add_argument("--ensembl_fasta", nargs="+", required=True)
        parser.add_argument("--ensembl_fasta_regex", default="(\w*)\.?", required=False)
        parser.add_argument("--refseq_fasta", nargs="+", default=[], required=False)
        parser.add_argument("--refseq_fasta_regex", default="(.*?)\s", required=False)
        parser.add_argument("--output", required=True)
        args = parser.parse_args()

        print('Getting uniprot FASTA sequences...')
        uniprot_swissprot = get_fasta_sequences(fasta_files=args.uniprot_fasta,
                                                    matching_regex=".*?\\|(.*?)\\|",
                                                    filter_function=lambda x: x.startswith("sp"),
                                                    reversed=True)
        uniprot_trembl = get_fasta_sequences(fasta_files=args.uniprot_fasta,
                                                    matching_regex=".*?\\|(.*?)\\|",
                                                    filter_function=lambda x: x.startswith("tr"),
                                                    reversed=True)

        print('TREMBL finished. \n Getting RefSeq sequences...')
        refseq_sequences = get_fasta_sequences(fasta_files=args.refseq_fasta,
                                                matching_regex=args.refseq_fasta_regex,
                                                reversed=True)
        print('RefSeq finished. \n Getting ENSEMBL sequences...')
        ensembl_proteins = get_fasta_sequences(fasta_files=args.ensembl_fasta,
                                                matching_regex=args.ensembl_fasta_regex,
                                                reversed=True)

        
        print('ENSEMBL finished. \n Merging...')
        all_proteins = dict(merge_fasta_dicts({
                "swissprot": uniprot_swissprot,
                "trembl": uniprot_trembl,
                "refseq": refseq_sequences,
                "ensembl": ensembl_proteins
            }))

        with open(args.output, "w") as f :
            json.dump(all_proteins, f)

    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        traceback.print_exc()
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)

if __name__ == "__main__":
    main()
