#!/usr/bin/env python3
"""
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez

Script that reformats Nonsense-Mediated Decay information
from SQANTI outputs to GTF.
"""
import argparse, sys, os, logging, csv
from IsoAnnot import openfile
from t2goAnnotationFile import get_structural_classification_df

def get_nmd_info(filename):
    '''
    Get Nonsense Mediated Decay information 
    '''
    with openfile(filename, "rt") as read_handler:
        tsv_reader = csv.reader(read_handler, delimiter="\t")
        for line_fields in tsv_reader:
            if line_fields[1] == "NMD":
                yield line_fields[0]

def main():
    try:
        parser = argparse.ArgumentParser(description='IsoAnnot NMD Layer')
        parser.add_argument('--nmd_file', nargs=None, required=True)
        parser.add_argument('--classification_file', nargs=None, required=True)
        parser.add_argument('--output', nargs=None, required=True)

        args = parser.parse_args()

        # Get info from files
        transcript_info = get_structural_classification_df(args.classification_file)
        nmd_info = get_nmd_info(args.nmd_file)

        with open(args.output, "w") as output_file:
            logging.info(f"Starting NMD conversion to GTF {args.output}")
            output_file_tsv = csv.writer(output_file, delimiter="\t")

            for transcript_id in nmd_info:
                transcript_params = transcript_info.get(transcript_id)
                output_file_tsv.writerow([
                        transcript_id,  # seqname
                        "NMD",  # source
                        "NMD",  # feature
                        ".",  # start
                        ".",  # end
                        ".",  # score
                        transcript_params.get("strand"),  # strand
                        ".",  # frame
                        f"ID=NMD; Name=nonsense-mediated decay"  # attribute
                    ])


    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
