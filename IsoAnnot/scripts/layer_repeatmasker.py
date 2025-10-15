#!/usr/bin/env python3

'''
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez
Modified by: Alessandra Martínez

Script that parses RepeatMasker's output and reformats
it into GTF.
'''

import argparse, sys, os, logging, csv, itertools
from IsoAnnot import openfile
from t2goAnnotationFile import get_structural_classification_df

def get_repeatmasker_info(filename):
    '''
    Retrieves RepeatMasker information
    '''
    with openfile(filename, "rt") as read_handler:
        tsv_reader = csv.reader(itertools.islice(read_handler, 3, None), delimiter=" ", skipinitialspace=True)
        for line_fields in tsv_reader:
            line_dict = {
                "identifier": line_fields[10],
                "transcript": line_fields[4],
                "start": int(line_fields[5]),
                "end": int(line_fields[6]),
                "desc": line_fields[9]
            }

            yield line_dict


def main():
    try:
        parser = argparse.ArgumentParser(description='IsoAnnot RepeatMasker Layer')
        parser.add_argument('--repeatmasker_file', nargs=None, required=True)
        parser.add_argument('--keep_version', nargs=None, default="False")
        parser.add_argument('--classification_file', nargs=None, required=True)
        parser.add_argument('--output', nargs=None, required=True)

        args = parser.parse_args()

        # Get information from files
        transcript_info = get_structural_classification_df(args.classification_file)
        repeatmasker_info = get_repeatmasker_info(args.repeatmasker_file)


        with open(args.output, "w") as output_file:
            logging.info(f"Starting RepeatMasker conversion to GTF {args.output}")
            output_file_tsv = csv.writer(output_file, delimiter="\t")

            for transcript_repeatmasker in repeatmasker_info:
                transcript_id = transcript_repeatmasker.get("transcript")

                if args.keep_version.upper()=="FALSE": 
                    transcript_id = transcript_id.split(".")[0]

                if transcript_id in transcript_info:
                    transcript_params = transcript_info.get(transcript_id)
                    repeat_id = transcript_repeatmasker.get("identifier")
                    repeat_desc = transcript_repeatmasker.get("desc")
                    repeat_start =  transcript_repeatmasker.get("start")
                    repeat_end = transcript_repeatmasker.get("end")

                    # Correct the positions
                    if repeat_start > repeat_end:
                        repeat_start, repeat_end = repeat_end, repeat_start

                    output_file_tsv.writerow([
                            transcript_id,  # seqname
                            "RepeatMasker",  # source
                            "repeat",  # feature
                            repeat_start,  # start
                            repeat_end,  # end
                            ".",  # score
                            transcript_params.get("strand"),  # strand
                            ".",  # frame
                            f"ID={repeat_id}; Name={repeat_id}; Desc={repeat_desc}"  # attribute
                        ])


    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
