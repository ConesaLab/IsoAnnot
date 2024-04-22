#!/usr/bin/env python3
'''
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez
Modified by: Alessandra Martínez

Script that retrieves Reactome information at gene resolution
(all the transcripts coming from the same gene will have the same
Reactome entries).
'''

import argparse, sys, os, logging, csv
import pandas as pd
from IsoAnnot import openfile
from t2goAnnotationFile import get_structural_classification_df
from pybiomart import Dataset


def get_reactome_info(filename):
    '''
    Retrieve information from reactome file.
    '''
    with openfile(filename, "rt") as read_handler:
        tsv_reader = csv.reader(read_handler, delimiter="\t")

        for line_fields in tsv_reader:
            line_dict = {
                "gene": line_fields[0],
                "id": line_fields[1],
                "name": line_fields[3]
            }
            yield line_dict


def write_output(transcript_id, reactome_id, pathway, output_filename):
    output_filename.writerow([
            transcript_id,  # seqname
            "REACTOME",  # source
            "PATHWAY",  # feature
            ".",  # start
            ".",  # end
            ".",  # score
            ".",  # strand
            ".",  # frame
            f"ID={reactome_id}; Name={pathway}"  # attribute
        ])


def main():
    try:
        parser = argparse.ArgumentParser(description='IsoAnnot Reactome Layer')
        parser.add_argument('--reactome_file', nargs=None, required=True)
        parser.add_argument('--classification_file', nargs=None, required=True)
        parser.add_argument('--output', nargs=None, required=True)
        parser.add_argument('--biomart_host', nargs="?", const="http://www.ensembl.org")
        parser.add_argument('--biomart_dataset', nargs=None, required=True)
        parser.add_argument('--species', nargs=None, required=True)

        args = parser.parse_args()

        if "plants" in args.biomart_host:
            biomart_schema = "plants_mart"
        else:
            biomart_schema = "default"

        # Get info from files
        transcript_info = get_structural_classification_df(args.classification_file)
        reactome_info = pd.read_csv(args.reactome_file, sep='\t', header=None)
        reactome_info.rename(columns={reactome_info.columns[0]: 'ensembl_gene_id', reactome_info.columns[1]: 'reactome_ID', reactome_info.columns[2]: 'URL', reactome_info.columns[3]: 'pathway', reactome_info.columns[4]: 'evidence', reactome_info.columns[5]: 'species'}, inplace=True)
        reactome_info = reactome_info[reactome_info.species == args.species]

        # Get ENSG-Gene_name conversion with Biomart and connect it with
        # the reactome dataframe
        dataset = Dataset(name=args.biomart_dataset, host=args.biomart_host, virtual_schema=biomart_schema)
        gene_name_table = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
        gene_name_table.rename(columns={'Gene stable ID': 'ensembl_gene_id', 'Gene name': 'gene_name'}, inplace=True)
        gene_name_table = gene_name_table.dropna()
        reactome_info = pd.merge(reactome_info, gene_name_table, how='inner', on='ensembl_gene_id')

        # For each transcript, annotate Reactome entries linked to its associated gene 
        with open(args.output, "w") as output_file:
            logging.info(f"Starting Reactome conversion to GTF {args.output}")
            output_file_tsv = csv.writer(output_file, delimiter="\t")
            for transcript, info in transcript_info.items():
                select_gene = reactome_info.loc[reactome_info["gene_name"] == info["ref_gene_name"]]
                if select_gene.empty:
                    continue
                select_gene_dict = select_gene.to_dict(orient="records")
                for record in select_gene_dict:  
                    reactome_id = record["reactome_ID"]
                    pathway = record["pathway"]
                    write_output(transcript,reactome_id, pathway, output_file_tsv)


    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
