#!/usr/bin/env python3

"""
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez

Script to retrieve information from Interproscan.
"""

import argparse, sys, os, logging, csv
from IsoAnnot import openfile

def read_interproscan_file(filename):
    """
    Retrieves information from Interproscan file

    Args:
        filename (str)
    
    Returns: 
        line_dict (dict)
    """
    dict_feature = {"PFAM": "DOMAIN",
                    "SIGNALP_EUK": "SIGNAL",
                    "MOBIDB_LITE": "DISORDER",
                    "TMHMM": "TRANSMEM",
                    "COILS": "COILED"}

    dict_desc = {"PFAM": "pfam_domain",
                 "SIGNALP_EUK": "signal peptide cleavage site",
                 "MOBIDB_LITE": "consensus disorder prediction",
                 "TMHMM": "Region of a membrane-bound protein predicted to be embeded in the membrane",
                 "COILS": "coiled-coil"}

    with openfile(filename, "rt") as read_handler:
        tsv_reader = csv.reader(read_handler, delimiter="\t")

        for line_fields in tsv_reader:
            line_dict = {
                "transcript": line_fields[0],
                "database": line_fields[1].upper(),
                "id": line_fields[3],
                "name": line_fields[4], 
                "feature": dict_feature.get(line_fields[1].upper()),
                "desc": line_fields[5] if line_fields[5] != "None" else dict_desc.get(line_fields[5]),
                "start": line_fields[6],
                "end": line_fields[7]
            }

            yield line_dict


def get_pfam_classification(filename):
    """
        Retrieves information from Pfam file

    Args:
        filename (str)
    
    Returns: 
        output_dict (dict)
    """
    output_dict = {}

    with openfile(filename, "rt") as read_handler:
        tsv_reader = csv.reader(read_handler, delimiter="\t")

        for line_fields in tsv_reader:
            pfam = line_fields[0]
            output_dict[pfam] = {
                "id": line_fields[1],
                "name": line_fields[3],
                "desc": line_fields[4]
            }

    return output_dict


def get_structural_classification_prot(classification_filename):
    '''
    Gets structural information for each transcript from the 
    SQANTI classification file.
    '''
    structural_classification = {}
    map_transcript_protein = {}

    with openfile(classification_filename, "rt") as read_handler:
        tsv_reader = csv.DictReader(read_handler, delimiter="\t")

        for line in tsv_reader:
            transcript_id = line.get("isoform")
            if line.get("coding") == "coding" and "refProt" in line:
                transcript_id = line.get("refProt")
            transcript_location = None

            transcript_qualifiers = {
                "length": line.get("length"),
                "chrom": line.get("chrom"),
                "strand": line.get("strand"),
                "ref_transcript": line.get("associated_transcript"),
                "ref_gene_name": line.get("associated_gene").replace(",", "_"),
                "primary_class": line.get("structural_category"),
                "second_class": ["RT_switching"] if line.get("RTS_stage") == "positive" else [],
                "cds_start": None if line.get("coding") != "coding" else line.get("CDS_start"),
                "cds_end": None if line.get("coding") != "coding" else line.get("CDS_end"),
                "ORF_length": None if line.get("coding") != "coding" else line.get("ORF_length"),
                "isoform": line.get("isoform"),
                "ref_protein": None,
                "ref_CDS": None,
                "extra_ref_protein": set([])
            }

            # Override CDS qualifiers for coding lines
            if line.get("coding") == "coding" and "refProt" in line:
                transcript_qualifiers["ref_protein"] = line.get("refProt")
                transcript_qualifiers["ref_CDS"] = line.get("refProt")
                map_transcript_protein[line.get("refProt")] = transcript_id

            structural_classification[transcript_id] = transcript_qualifiers

    return structural_classification


def main():
    try:
        parser = argparse.ArgumentParser(description='IsoAnnot InterProScan Layer')
        parser.add_argument('--interproscan_file', nargs=None, required=True)
        parser.add_argument('--pfam_file', nargs=None, required=True)
        parser.add_argument('--keep_version', nargs=None, default="False")
        parser.add_argument('--classification_file', nargs=None, required=True)
        parser.add_argument('--allowed_databases', nargs=1, default=["PFAM", "SIGNALP_EUK", "MOBIDB_LITE", "COILS", "TMHMM"])
        parser.add_argument('--output', nargs=None, required=True)

        args = parser.parse_args()

        # Get info from files
        pfam_classification = get_pfam_classification(args.pfam_file)
        transcript_info = get_structural_classification_prot(args.classification_file)
        interproscan_info = read_interproscan_file(args.interproscan_file)

        with open(args.output, "w") as output_file:
            logging.info(f"Starting InterProScan conversion to GTF {args.output}")
            output_file_tsv = csv.writer(output_file, delimiter="\t")
            
            # Get interproscan info
            for transcript_interproscan in interproscan_info:
                transcript_id = transcript_interproscan.get("transcript")
                transcript_db = transcript_interproscan.get("database")
                transcript_its_id = transcript_interproscan.get("id")
                
                # Keep transcript version or not
                if args.keep_version.upper() == "FALSE":
                    transcript_id = transcript_id.split(".")[0]
                
                # Only annotate features from the allowed databases
                if transcript_id in transcript_info and transcript_db in args.allowed_databases:
                    transcript_id = transcript_info.get(transcript_id).get("isoform")
                    transcript_name = transcript_interproscan.get("name")
                    transcript_desc = transcript_interproscan.get("desc")

                    output_file_tsv.writerow([
                        transcript_id,  # seqname
                        transcript_db,  # source
                        transcript_interproscan.get("feature"),  # feature
                        transcript_interproscan.get("start"),  # start
                        transcript_interproscan.get("end"),  # end
                        ".",  # score
                        ".",  # strand
                        ".",  # frame
                        f"ID={transcript_its_id}; Name={transcript_name}; Desc={transcript_desc}" # attribute
                    ])

                    # Add additional rows for PFAM database
                    if transcript_db == "PFAM" and transcript_its_id in pfam_classification:
                        pfam_id = pfam_classification.get(transcript_its_id).get("id")
                        pfam_name = pfam_classification.get(transcript_its_id).get("name")
                        pfam_desc = pfam_classification.get(transcript_its_id).get("desc")

                        output_file_tsv.writerow([
                            transcript_id,  # seqname
                            transcript_db,  # source
                            "clan",  # feature
                            ".",  # start
                            ".",  # end
                            ".",  # score
                            ".",  # strand
                            ".",  # frame
                            f"ID={pfam_id}; Name={pfam_name}; Desc={pfam_desc}"  # attribute
                        ])

    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
