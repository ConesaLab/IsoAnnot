#!/usr/bin/env python3
'''
Script that links each coding transcript of the reference (ensembl and refseq modes) to every protein
ID with the same sequence. The output has the protein_assoc format of transcript2reference.py
(mytranscripts), and t2goAnnotationFile.py uses it to choose the protein ID of the GFF3 (GitHub #10).
'''

import argparse, csv, json, logging, os, sys, traceback
from IsoAnnot import openfile

# The legacy tappAS files only list UniProt and the native database of the mode.
ALLOWED_DATABASES = {"ensembl": ["swissprot", "trembl", "ensembl"],
                     "refseq": ["swissprot", "trembl", "refseq"]}


def main():
    try:
        parser = argparse.ArgumentParser(description="Reference protein association")
        parser.add_argument("--classification_file", required=True)
        parser.add_argument("--species_db", required=True, help="parsed_databases.json from parse_protein_databases.py")
        parser.add_argument("--database", required=True, choices=sorted(ALLOWED_DATABASES))
        parser.add_argument("--output", required=True)
        args = parser.parse_args()

        with open(args.species_db) as fp:
            all_proteins = json.load(fp)

        # parsed_databases.json maps sequence -> IDs. Invert it for the native IDs. An ID that maps to
        # several sequences once its version is removed (AT1G01010.1/.2 in Ensembl Plants) can't be used.
        native_sequences = {}
        for seq, protein_dbs in all_proteins.items():
            for protein_id in protein_dbs.get(args.database, []):
                native_sequences[protein_id] = seq if native_sequences.get(protein_id, seq) == seq else None

        n_coding = 0
        n_missing = 0
        with openfile(args.classification_file, "rt") as read_handler, open(args.output + ".tmp", "w") as output_file:
            output_tsv = csv.writer(output_file, delimiter="\t", lineterminator="\n")
            for line in csv.DictReader(read_handler, delimiter="\t"):
                if line["coding"] != "coding" or line["refProt"] == "NA":
                    continue
                n_coding += 1

                # The json has unversioned IDs (NP_000005.3 -> NP_000005), as get_fasta_sequences writes them
                seq = native_sequences.get(line["refProt"].split(".")[0])
                if seq is None:
                    n_missing += 1
                    continue

                for database, protein_ids in all_proteins[seq].items():
                    if database in ALLOWED_DATABASES[args.database]:
                        for protein_id in protein_ids:
                            output_tsv.writerow([line["isoform"], protein_id, database, "sequence"])

        os.replace(args.output + ".tmp", args.output)

        if n_missing:
            logging.warning(f"{n_missing} of {n_coding} coding transcripts have a refProt with no unique sequence in "
                            f"{args.species_db}. t2goAnnotationFile.py uses the refProt as their protein ID.")

    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        traceback.print_exc()
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
