#!/usr/bin/env python3

"""
Author: Carlos Martínez
Modified by: Alessandra Martinez

Script to retrieve Gene Ontology information at gene resolution
(all the transcripts coming from the same gene will have the same
GO terms).
"""

import argparse, sys, os, logging, csv, traceback
import pandas as pd
from pybiomart import Dataset
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
        f"ID={go_accession}; Name={go_name};"  # attribute
    ])


def main():

    allowed_features = ["molecular_function", "cellular_component", "biological_process"]
    try:
        parser = argparse.ArgumentParser(description='LayerGO')
        parser.add_argument('--classification_file', required=True)
        parser.add_argument('--biomart_host', nargs="?", const="http://www.ensembl.org")
        parser.add_argument('--biomart_dataset', nargs=None, required=True)
        parser.add_argument('--output', nargs=None, required=True)

        args = parser.parse_args()

        if "plants" in args.biomart_host:
            biomart_schema = "plants_mart"
        else:
            biomart_schema = "default"

        # Get info from files
        transcript_info = get_structural_classification_df(args.classification_file)

        # Open the connection to Biomart server and send query
        dataset = Dataset(name=args.biomart_dataset, host=args.biomart_host, 
                          virtual_schema=biomart_schema)
        go_table = dataset.query(attributes=["external_gene_name",
                                                "go_id",
                                                "name_1006", #go name
                                                "namespace_1003"]) #go category
        # Remove NAs and allow only specific go categories 
        go_table = go_table.dropna()
        go_table = go_table[go_table['GO domain'].isin(allowed_features)] 

        # For each transcript, annotate GO terms linked to its associated gene 
        with open(args.output, "w") as output_file:
            logging.info(f"Starting GO conversion to GTF {args.output}")
            output_file_tsv = csv.writer(output_file, delimiter="\t")
            for transcript, info in transcript_info.items():
                select_gene = go_table.loc[go_table["Gene name"] == info["ref_gene_name"]]
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
