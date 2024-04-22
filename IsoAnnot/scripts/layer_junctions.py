#!/usr/bin/env python3
"""
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez

Script that reformats junctions information from SQANTI output.
"""

import argparse, sys, os, logging, csv
from IsoAnnot import openfile
from t2goAnnotationFile import get_structural_classification_df


def get_junctions_info(filename):
    """
    Get junctions info from file line by line

    Args:
        filename (str)

    Returns:
        line (list): junctions line
    """
    with openfile(filename, "rt") as read_handler:
        tsv_reader = csv.DictReader(read_handler, delimiter="\t")
        for line in tsv_reader:
            yield line


def main():
    try:
        parser = argparse.ArgumentParser(description='IsoAnnot Junctions Layer')
        parser.add_argument('--junctions_file', nargs=None, required=True)
        parser.add_argument('--classification_file', nargs=None, required=True)
        parser.add_argument('--output', nargs=None, required=True)

        args = parser.parse_args()

        # Get info from files
        transcript_info = get_structural_classification_df(args.classification_file)
        junctions_info = get_junctions_info(args.junctions_file)

        with open(args.output, "w") as output_file:
            logging.info(f"Starting junctions conversion to GTF {args.output}")
            output_file_tsv = csv.writer(output_file, delimiter="\t")

            for transcript_junction in junctions_info:
                transcript_id = transcript_junction.get("isoform")
                transcript_params = transcript_info.get(transcript_id)

                gtf_attributes = f"ID={transcript_junction['junction_category']}_{transcript_junction['canonical']}; Chr=chr{transcript_params['chrom']}"

                output_file_tsv.writerow([
                        transcript_id,  # seqname
                        "tappAS",  # source
                        "splice_junction",  # feature
                        transcript_junction.get("genomic_start_coord"),  # start
                        transcript_junction.get("genomic_end_coord"),  # end
                        ".",  # score
                        transcript_info.get(transcript_id).get("strand"),  # strand
                        ".",  # frame
                        gtf_attributes  # attribute
                    ])

    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
