#!/usr/bin/env python3
'''
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez 
Modified by: Alessandra Martinez

Script that parses UTRscan output and writes it in
GTF format
'''

import argparse, sys, os, re, logging, csv, traceback
from IsoAnnot import openfile
from t2goAnnotationFile import get_structural_classification_df


def read_utrscan_files(filename):
    motif_regex =  re.compile(r"^(U[0-9]+)\s+-\s([\w-]+)\s+-\s+(.+)$")
    transcript_regex = re.compile(r"^(.+?)\s.*:\s([\w-]+)\s\[([0-9]+),([0-9]+)\]")

    motif_info = {}

    with openfile(filename, "rt") as read_handler:
        for line in read_handler:
            # First lines of UTRscan output contain the definition of domain.
            # Extract the info from them.
            if line.startswith("U"):
                motif_capture = re.match(motif_regex, line)
                motif_id, motif, motif_desc = motif_capture.groups()
                motif_info[motif] = {"id": motif_id, "desc": motif_desc}

            # The next lines contain the results of UTRScan for each transcript
            # TODO: is enough with this check?
            elif ":" in line:
                transcript_info = re.match(transcript_regex, line)
                transcript, motif, cord_start, cord_end = transcript_info.groups()

                line_dict = {
                    "transcript": transcript,
                    "motif": motif,
                    "motif_id": motif_info.get(motif).get('id'),
                    "motif_desc": motif_info.get(motif).get('desc'),
                    "start": cord_start,
                    "end": cord_end
                }
                yield line_dict


def main():

    try: 

        parser = argparse.ArgumentParser(description='UTRScan')
        parser.add_argument('--utrscan_file', nargs=None, required=True)
        parser.add_argument('--keep_version', nargs=None, default="False")
        parser.add_argument('--classification_file', nargs=None, required=True)
        parser.add_argument('--output', nargs=None, required=True)

        args = parser.parse_args()

        # Get info from files
        transcript_info = get_structural_classification_df(args.classification_file)
        utrscan_info = read_utrscan_files(args.utrscan_file)

        with open(args.output, "w") as output_file:
            logging.info(f"Starting UTRScan conversion to GTF {args.output}")

            output_file_tsv = csv.writer(output_file, delimiter="\t")
            for transcript_utrscan in utrscan_info:
                transcript_id = transcript_utrscan.get("transcript")

                if args.keep_version.upper()=="FALSE":
                    transcript_id = transcript_id.partition(".")[0]

                if transcript_id in transcript_info:
                    transcript_params = transcript_info.get(transcript_id)

                    # 5'UTR
                    if int(transcript_utrscan.get("end")) < int(transcript_params.get("cds_start") or 0):
                        gtf_attributes = f"ID={transcript_utrscan['motif_id']}; Name={transcript_utrscan['motif']}; Desc={transcript_utrscan['motif_desc']}"
                        output_file_tsv.writerow([
                            transcript_id,  # seqname
                            "UTRsite",  # source
                            "uORF" if transcript_utrscan.get("motif") == "uORF" else "5UTRmotif",
                            transcript_utrscan.get("start"),  # start
                            transcript_utrscan.get("end"),  # end
                            ".",  # score
                            transcript_params.get("strand"),  # strand
                            ".",  # frame
                            gtf_attributes  # attribute
                        ])
                        
                    # 3'UTR
                    elif int(transcript_params.get("cds_end") or 0) < int(transcript_utrscan.get("start")):
                        gtf_attributes = f"ID={transcript_utrscan['motif_id']}; Name={transcript_utrscan['motif']}; Desc={transcript_utrscan['motif_desc']}"
                        output_file_tsv.writerow([
                            transcript_id,  # seqname
                            "UTRsite",  # source
                            "PAS" if "PAS" in transcript_utrscan.get("motif") else "3UTRmotif",
                            transcript_utrscan.get("start"),  # start
                            transcript_utrscan.get("end"),  # end
                            ".",  # score
                            transcript_params.get("strand"),  # strand
                            ".",  # frame
                            gtf_attributes  # attribute
                        ])

    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        traceback.print_exc()
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)

if __name__ == "__main__":
    main()
