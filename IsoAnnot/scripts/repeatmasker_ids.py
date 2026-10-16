#!/usr/bin/env python3
"""
RepeatMasker stops when a FASTA ID is longer than 50 characters, which long-read isoform IDs can be.
"shorten" writes a copy of the FASTA where those IDs are replaced by their MD5 hash, plus a hash -> ID table.
"restore" puts the original IDs back in the RepeatMasker .out file. Used by rule run_repeatmasker.
"""

import argparse, hashlib, logging, os
from IsoAnnot import openfile

# Longest sequence ID RepeatMasker accepts
MAX_ID_LENGTH = 50


def shorten_fasta_ids(fasta, out_fasta, mapping_file):
    """
    Copy fasta to out_fasta, replacing IDs longer than MAX_ID_LENGTH by their MD5 hash (32 characters).
    Headers of shortened IDs lose their description, which RepeatMasker does not use.

    :return: number of IDs replaced. mapping_file (hash, original ID) is always written, empty if 0.
    """
    n_shortened = 0
    with openfile(fasta, "rt") as f_in, open(f"{out_fasta}.tmp", "w") as f_out, open(f"{mapping_file}.tmp", "w") as mapping:
        for line in f_in:
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                if len(seq_id) > MAX_ID_LENGTH:
                    short_id = hashlib.md5(seq_id.encode()).hexdigest()
                    mapping.write(f"{short_id}\t{seq_id}\n")
                    line = f">{short_id}\n"
                    n_shortened += 1
            f_out.write(line)
    os.replace(f"{out_fasta}.tmp", out_fasta)
    os.replace(f"{mapping_file}.tmp", mapping_file)
    return n_shortened


def restore_out_ids(repeatmasker_file, mapping_file, output):
    """
    Write the RepeatMasker .out file to output with the hashed IDs of the query sequence column (5th)
    replaced by the original IDs. Header lines and the "no repetitive sequences" line are copied as they are.

    :return: number of hit lines whose ID was restored.
    """
    with open(mapping_file) as mapping:
        original_ids = dict(line.rstrip("\n").split("\t") for line in mapping)

    n_restored = 0
    with open(repeatmasker_file) as f_in, open(f"{output}.tmp", "w") as f_out:
        for line in f_in:
            fields = line.split()
            if len(fields) > 4 and fields[4] in original_ids:
                # The four columns before the query are numbers, so the first match is the query column
                line = line.replace(fields[4], original_ids[fields[4]], 1)
                n_restored += 1
            f_out.write(line)
    os.replace(f"{output}.tmp", output)
    return n_restored


def main():
    parser = argparse.ArgumentParser(description="Shorten long FASTA IDs for RepeatMasker and restore them in its .out file")
    subparsers = parser.add_subparsers(dest="command", required=True)

    shorten = subparsers.add_parser("shorten", help="Write a FASTA with IDs longer than 50 characters hashed")
    shorten.add_argument("--fasta", required=True)
    shorten.add_argument("--output", required=True, help="FASTA to give to RepeatMasker")
    shorten.add_argument("--mapping", required=True, help="Hash -> original ID table")

    restore = subparsers.add_parser("restore", help="Put the original IDs back in a RepeatMasker .out file")
    restore.add_argument("--repeatmasker_file", required=True)
    restore.add_argument("--mapping", required=True, help="Table written by shorten")
    restore.add_argument("--output", required=True)

    args = parser.parse_args()

    if args.command == "shorten":
        n_shortened = shorten_fasta_ids(args.fasta, args.output, args.mapping)
        logging.info(f"{n_shortened} IDs longer than {MAX_ID_LENGTH} characters were replaced by their MD5 hash. They are listed in {args.mapping}")
    else:
        n_restored = restore_out_ids(args.repeatmasker_file, args.mapping, args.output)
        logging.info(f"Original IDs restored in {n_restored} RepeatMasker hits")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
