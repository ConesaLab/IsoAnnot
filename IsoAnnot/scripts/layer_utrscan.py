#!/usr/bin/env python3
'''
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez 
Modified by: Alessandra Martinez

Script that parses UTRscan output and writes it in
GTF format
'''

import argparse, sys, os, re, logging, csv, traceback
from collections import Counter
from IsoAnnot import openfile, normalize_id
from t2goAnnotationFile import get_structural_classification_df


UTRSITE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "software", "bin", "UtrSite.Command")


def read_utrsite_regions(filename):
    """
    UTR where each UTRsite motif acts, from the catalogue UtrScan runs with: "5", "3" or "all".
    Its headers are "// Signal : uORF|U0033|Upstream Open Reading Frame (uORF)|5|Pattern Similarity|*2759*".
    """
    regions = {}
    with open(filename) as read_handler:
        for line in read_handler:
            if line.startswith("// Signal"):
                fields = line.split(":", 1)[1].split("|")
                regions[fields[0].strip()] = fields[3].strip()
    return regions


def get_motif_feature(motif, motif_region, start, end, cds_start, cds_end, keep_noncoding):
    """
    Feature of a UtrScan hit, or None when it is not written. UtrScan searches every motif along the whole
    transcript, so a motif is kept only in the UTR its catalogue entry gives: an "upstream ORF" after the
    CDS is not one. Transcripts without a CDS have no UTR: their motifs are written as ncRNAmotif, without
    the ones of the 5'UTR (uORF, IRES...), which act on the translation of a CDS.

    :param motif_region: "5", "3" or "all". Motifs that are not in the catalogue count as "all".
    :param cds_start, cds_end: None for a non-coding transcript
    """
    if cds_start is None:
        if not keep_noncoding or motif_region == "5":
            return None
        return "PAS" if "PAS" in motif else "ncRNAmotif"
    if end < cds_start:
        if motif_region == "3":
            return None
        return "uORF" if motif == "uORF" else "5UTRmotif"
    if cds_end < start:
        if motif_region == "5":
            return None
        return "PAS" if "PAS" in motif else "3UTRmotif"
    return None


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
        parser.add_argument('--classification_file', nargs=None, required=True)
        parser.add_argument('--output', nargs=None, required=True)
        parser.add_argument('--utrsite_file', default=UTRSITE_FILE, help="UTRsite catalogue UtrScan ran with (default: %(default)s)")
        parser.add_argument('--skip_noncoding', action="store_true", help="Leave out the motifs of transcripts without a CDS")

        args = parser.parse_args()

        # Get info from files
        transcript_info = get_structural_classification_df(args.classification_file)
        utrscan_info = read_utrscan_files(args.utrscan_file)
        motif_regions = read_utrsite_regions(args.utrsite_file)
        n_features = Counter()

        with open(args.output, "w") as output_file:
            logging.info(f"Starting UTRScan conversion to GTF {args.output}")

            output_file_tsv = csv.writer(output_file, delimiter="\t")
            for transcript_utrscan in utrscan_info:
                transcript_id = transcript_utrscan.get("transcript")

                # Strip the version only when the exact ID is not in the classification. transcript_versioned is
                # not used here, so a wrong value in the config cannot empty the layer.
                if transcript_id not in transcript_info:
                    transcript_id = normalize_id(transcript_id)

                if transcript_id in transcript_info:
                    transcript_params = transcript_info.get(transcript_id)
                    coding = transcript_params.get("cds_start") is not None
                    feature = get_motif_feature(transcript_utrscan["motif"], motif_regions.get(transcript_utrscan["motif"], "all"),
                                                int(transcript_utrscan["start"]), int(transcript_utrscan["end"]),
                                                int(transcript_params["cds_start"]) if coding else None,
                                                int(transcript_params["cds_end"]) if coding else None,
                                                keep_noncoding=not args.skip_noncoding)
                    n_features[feature or ("left out, coding" if coding else "left out, non-coding")] += 1
                    if feature:
                        gtf_attributes = f"ID={transcript_utrscan['motif_id']}; Name={transcript_utrscan['motif']}; Desc={transcript_utrscan['motif_desc']}"
                        output_file_tsv.writerow([
                            transcript_id,  # seqname
                            "UTRsite",  # source
                            feature,
                            transcript_utrscan.get("start"),  # start
                            transcript_utrscan.get("end"),  # end
                            ".",  # score
                            transcript_params.get("strand"),  # strand
                            ".",  # frame
                            gtf_attributes  # attribute
                        ])

        logging.info("UtrScan hits: " + ", ".join(f"{count} {name}" for name, count in sorted(n_features.items())))

    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        traceback.print_exc()
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
