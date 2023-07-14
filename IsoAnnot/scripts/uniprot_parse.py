#!/usr/bin/env python3
'''
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez
Modified by: Alessandra Martinez

Script that takes Complete SwissProt and UniProtKB/TrEMBL data set  in flat
file format and extracts relevant information.
'''

from Bio import SwissProt
import argparse, sys, os, csv, logging, traceback
from IsoAnnot import openfile


def main():
    parser = argparse.ArgumentParser(description="UniprotParser")
    parser.add_argument("--uniprot_files", nargs="+", required=True) 
    parser.add_argument("--output", required=True, type=str)

    args = parser.parse_args()

    cross_ref_information = ("RefSeq", "Ensembl", "EnsemblPlants", "EnsemblMetazoa", )

    try:
        with open(args.output, "w", newline="") as output_handler: 
            tsv_writer = csv.writer(output_handler, delimiter='\t')
            
            # Parse both SwissProt and TrEMBL datasets in flat file format
            for input_file in args.uniprot_files:  
                logging.info(f"Starting analysis of {input_file}...\n")

                with openfile(input_file, "rb") as read_handler: 
                    for record in SwissProt.parse(read_handler):  
                        taxa = record.taxonomy_id[0]
                        accession = record.accessions[0]
                        cross_references = record.cross_references
                        record_features = record.features

                        # Write Ensembl and Refseq cross-references
                        for reference in cross_references:
                            if reference[0] in cross_ref_information:
                                tsv_writer.writerow([accession, taxa, "DR", *reference])

                        # Write features
                        for feature in record_features:                            
                            ftype = feature.type
                            fstart = str(feature.location).split(":")[0][1:]
                            # FT only present in specific isoforms
                            # TODO they are not parsed correctly
                            if "[" in fstart:
                                continue
                            fend = str(feature.location).split(":")[1][0:-1]  
                            if fstart == "UnknownPosition()" or fend == "UnknownPosition()":
                                continue
                            lfqualifiers = []
                            for i in feature.qualifiers:
                                lfqualifiers.append(feature.qualifiers[i])
                            fqualifiers = ". ".join(lfqualifiers)
                            tsv_writer.writerow([accession, taxa, "FT", ftype, fstart, fend, fqualifiers]) 


                logging.info(f"Finished analysis of {input_file}...\n")
    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        traceback.print_exc()
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)

if __name__ == "__main__":
    main()

