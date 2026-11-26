#!/usr/bin/env python3

"""
Author: Carlos Martínez
Modified by: Alessandra Martinez
Modified by: Raúl López

Script to retrieve Gene Ontology information at gene resolution
(all the transcripts coming from the same gene will have the same
GO terms).
"""

import argparse, sys, os, logging, csv, traceback
import pandas as pd
from IsoAnnot import read_biomart_table
from t2goAnnotationFile import get_structural_classification_df


def write_output(transcript_id, go_feature, go_accession, go_name, output_filename):
    output_filename.writerow([
        transcript_id,  # seqname
        "GeneOntology",  # source
        go_feature.split("_")[1][0].upper(),  # feature            
        ".",  # start
        ".",  # end
        ".",  # score
        ".",  # strand
        ".",  # frame
        f"ID={go_accession}; Name={go_name}"  # attribute
    ])


def main():

    allowed_features = ["molecular_function", "cellular_component", "biological_process"]
    try:
        parser = argparse.ArgumentParser(description='LayerGO')
        parser.add_argument('--classification_file', required=True)
        parser.add_argument('--output', nargs=None, required=True)
        parser.add_argument('--db', nargs=None, required=True)
        parser.add_argument('--biomart_table', required=True, help='Table written by save_biomart_tables.py, with the BioMart answer of the Ensembl release in use')
        args = parser.parse_args()

        # Get info from files
        transcript_info = get_structural_classification_df(args.classification_file)

        # Select the appropriate gene identifier column based on the database choice
        gene_column = "Gene stable ID" if args.db == "mytranscripts" else "Gene name"

        # The table has one row per gene ID, name and GO term: the rows of the columns used here are made unique
        go_columns = [gene_column, "GO term accession", "GO term name", "GO domain"]
        go_table = read_biomart_table(args.biomart_table, "layer_go")[go_columns].drop_duplicates()
        # Remove NAs and allow only specific go categories 
        go_table = go_table.dropna()
        go_table = go_table[go_table['GO domain'].isin(allowed_features)] 

        # For each transcript, annotate GO terms linked to its associated gene 
        with open(args.output, "w") as output_file:
            logging.info(f"Starting GO conversion to GTF {args.output}")
            output_file_tsv = csv.writer(output_file, delimiter="\t")
            for transcript, info in transcript_info.items():
                select_gene = go_table.loc[go_table[gene_column] == info["ref_gene_name"]]
                if select_gene.empty:
                    continue
                select_gene_dict = select_gene.to_dict(orient="records")
                for record in select_gene_dict:  
                    go_feature = record["GO domain"]
                    go_accession = record["GO term accession"]
                    go_name = record["GO term name"]
                    write_output(transcript, go_feature, go_accession, go_name, output_file_tsv)


    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        traceback.print_exc()
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()

