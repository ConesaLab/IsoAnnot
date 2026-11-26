#!/usr/bin/env python3
# Author: Lorena de la Fuente Lorente
# Refactored: Carlos Martínez
# Modified by Fran Pardo-Palacios and Alessandra Martinez

import argparse, csv, logging, sys, os, traceback, json
from IsoAnnot import merge_fasta_dicts, read_gtf_protein_ids, annotation_id_forms
from uniprotPhosphosite_genomicCoordinates import get_fasta_sequences


def get_native_sequences(fasta_files, matching_regex, gtf_file, database):
    """
    {sequence: [protein IDs]} of the Ensembl or RefSeq protein FASTA. The IDs are read whole and stored in
    the form the annotation of the database uses (annotation_id_forms), because these are the IDs
    written in the GFF3. Scripts that look them up remove the version on their side.
    """
    sequences = get_fasta_sequences(fasta_files=fasta_files, matching_regex=matching_regex,
                                    keep_version="True", reversed=True)
    if not gtf_file:
        return sequences
    annotation_ids = read_gtf_protein_ids(gtf_file)
    forms = annotation_id_forms([protein_id for protein_ids in sequences.values() for protein_id in protein_ids],
                                annotation_ids, database)
    return {sequence: [forms[protein_id] for protein_id in protein_ids] for sequence, protein_ids in sequences.items()}


def main():
    try: 
        parser = argparse.ArgumentParser(description="Parse Databases")
        parser.add_argument("--uniprot_fasta", nargs="+", required=True)
        parser.add_argument("--ensembl_fasta", nargs="+", required=True)
        # The whole first word: (\w*) stopped at the dot and read AT1G01010.1 (Ensembl Plants) as its gene
        parser.add_argument("--ensembl_fasta_regex", default=r"(\S+)", required=False)
        parser.add_argument("--ensembl_gtf", required=False, help="Ensembl annotation, for the form of its protein IDs")
        parser.add_argument("--refseq_fasta", nargs="+", default=[], required=False)
        parser.add_argument("--refseq_fasta_regex", default=r"(.*?)\s", required=False)
        parser.add_argument("--refseq_gtf", required=False, help="RefSeq annotation, for the form of its protein IDs")
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
        refseq_sequences = get_native_sequences(args.refseq_fasta, args.refseq_fasta_regex, args.refseq_gtf, "refseq")
        print('RefSeq finished. \n Getting ENSEMBL sequences...')
        ensembl_proteins = get_native_sequences(args.ensembl_fasta, args.ensembl_fasta_regex, args.ensembl_gtf, "ensembl")

        
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
    logging.basicConfig(level=logging.INFO)
    main()
