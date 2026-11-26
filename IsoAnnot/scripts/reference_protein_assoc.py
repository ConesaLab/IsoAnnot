#!/usr/bin/env python3
'''
Script that links each coding transcript of the reference (ensembl and refseq modes) to every protein
ID with the same sequence. The output has the protein_assoc format of transcript2reference.py
(mytranscripts), and t2goAnnotationFile.py uses it to choose the protein ID of the GFF3 (GitHub #10).
'''

import argparse, csv, json, logging, os, sys, traceback
from IsoAnnot import openfile, TRAILING_NUMBER_REGEX

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

        # parsed_databases.json maps sequence -> IDs, written as the annotation writes them. Invert it for
        # the native IDs.
        native_sequences = {}
        for seq, protein_dbs in all_proteins.items():
            for protein_id in protein_dbs.get(args.database, []):
                native_sequences[protein_id] = seq if native_sequences.get(protein_id, seq) == seq else None

        # The annotation and the protein FASTA can have different versions of a protein (NP_000005.3 and
        # NP_000005.4). Such a refProt is looked up without the final number, if that points to one sequence.
        sequences_without_number = {}
        for protein_id, seq in native_sequences.items():
            short_id = TRAILING_NUMBER_REGEX.sub("", protein_id)
            sequences_without_number[short_id] = seq if sequences_without_number.get(short_id, seq) == seq else None

        n_coding = 0
        n_missing = 0
        n_other_version = 0
        with openfile(args.classification_file, "rt") as read_handler, open(args.output + ".tmp", "w") as output_file:
            output_tsv = csv.writer(output_file, delimiter="\t", lineterminator="\n")
            for line in csv.DictReader(read_handler, delimiter="\t"):
                if line["coding"] != "coding" or line["refProt"] == "NA":
                    continue
                n_coding += 1

                seq = native_sequences.get(line["refProt"])
                if seq is None:
                    seq = sequences_without_number.get(TRAILING_NUMBER_REGEX.sub("", line["refProt"]))
                    n_other_version += seq is not None
                if seq is None:
                    n_missing += 1
                    continue

                for database, protein_ids in all_proteins[seq].items():
                    if database in ALLOWED_DATABASES[args.database]:
                        for protein_id in protein_ids:
                            output_tsv.writerow([line["isoform"], protein_id, database, "sequence"])

        os.replace(args.output + ".tmp", args.output)

        if n_other_version:
            logging.warning(f"{n_other_version} of {n_coding} coding transcripts have a refProt that {args.species_db} "
                            f"has with another version. Their protein IDs are the ones of the protein FASTA. Check "
                            f"that the annotation and the protein FASTA are of the same release.")
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
