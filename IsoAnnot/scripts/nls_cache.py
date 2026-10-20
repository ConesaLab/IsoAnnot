#!/usr/bin/env python3
"""
Store of NucImport results, so a protein sequence predicted in an earlier run is not predicted again.

NucImport scores each protein on its own and, with ID=F, from its sequence only: the same sequence gives
the same result whatever its ID or the other proteins in its chunk. A result is therefore stored under the
SHA-256 of the sequence, together with the settings it was computed with (model, organism, ID flag and
the MD5 of the jar), and only reused when all of them match.

The store is a folder of TSV files, one per run, holding only the sequences that run predicted. Parallel
runs then never write to the same file.

  lookup: splits the unique proteins of a run into the ones to predict and the ones already in the store.
  store:  adds the new NucImport results to the store, and stops the run if most proteins got no result.
"""

import argparse
import glob
import hashlib
import logging
import os
import sys
import time

from Bio import SeqIO

# First two lines of every NucImport output, as NucImportMay2012.jar writes them
NUCIMPORT_HEADER = ("Protein ID\tImport Probability\t cNLS probability\tPredicted Most probable cNLS\t"
                    "Predicted most probable cNLS class number\tPosition")
NUCIMPORT_STARS = "*" * 93

STORE_COLUMNS = ["seq_sha256", "model", "organism", "id_flag", "jar_md5",
                 "imp_prob", "cnls_prob", "cnls_seq", "cnls_class", "position"]

# Above this fraction of proteins without a result, NucImport did not work and the run stops
MAX_MISSING_FRACTION = 0.5


def sequence_hash(seq):
    return hashlib.sha256(seq.encode()).hexdigest()


def file_md5(path):
    md5 = hashlib.md5()
    with open(path, "rb") as h:
        for block in iter(lambda: h.read(1 << 20), b""):
            md5.update(block)
    return md5.hexdigest()


def read_store(cache_dir, settings):
    """
    Reads every store file in cache_dir and keeps the results computed with the given settings.

    :param settings: (model, organism, id_flag, jar_md5)
    :return: dict sequence hash -> list of the 5 NucImport values, as text
    """
    results = {}
    n_conflicts = 0
    for store_file in sorted(glob.glob(os.path.join(cache_dir, "*.tsv"))):
        with open(store_file) as h:
            header = h.readline().rstrip("\n").split("\t")
            if header != STORE_COLUMNS:
                logging.warning(f"{store_file} is not a NucImport store file and was skipped")
                continue
            for line in h:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != len(STORE_COLUMNS) or tuple(fields[1:5]) != settings:
                    continue
                values = fields[5:]
                # The first file wins. NucImport is deterministic, so two different values mean a problem
                if fields[0] in results:
                    if results[fields[0]] != values:
                        n_conflicts += 1
                    continue
                results[fields[0]] = values
    if n_conflicts:
        logging.warning(f"{n_conflicts} sequences have different results in different store files in {cache_dir}. "
                        "The first file in name order was used")
    return results


def read_nucimport_outputs(output_files):
    """
    :return: dict protein ID -> list of the 5 NucImport values, as text
    """
    results = {}
    for output_file in output_files:
        with open(output_file) as h:
            for line in h:
                if line.startswith("Protein ID\t") or line.startswith("*") or not line.strip():
                    continue
                # Every line ends with a tab, so it splits into 7 fields
                fields = line.rstrip("\n").split("\t")
                results[fields[0]] = fields[1:6]
    return results


def write_atomic(path, lines):
    """Writes to <path>.tmp and renames it, so an interrupted run never leaves a truncated file."""
    with open(f"{path}.tmp", "w") as out:
        out.writelines(lines)
    os.replace(f"{path}.tmp", path)


def lookup(args, settings):
    stored = {} if args.no_cache else read_store(args.cache_dir, settings)
    to_predict = []
    cached_lines = [f"{NUCIMPORT_HEADER}\n", f"{NUCIMPORT_STARS}\n"]
    n_proteins = 0
    for record in SeqIO.parse(args.fasta, "fasta"):
        n_proteins += 1
        values = stored.get(sequence_hash(str(record.seq)))
        if values is None:
            to_predict.append(f">{record.id}\n{record.seq}\n")
        else:
            # Same format as NucImport, under the protein ID of this run
            cached_lines.append("\t".join([record.id] + values) + "\t\n")

    write_atomic(args.to_predict, to_predict)
    write_atomic(args.cached, cached_lines)
    if args.no_cache:
        logging.info(f"NucImport cache turned off (nls_cache): all {n_proteins} proteins will be predicted")
    else:
        logging.info(f"{len(stored)} results in the NucImport cache {args.cache_dir}")
        logging.info(f"{n_proteins - len(to_predict)} of {n_proteins} proteins found in the cache, "
                     f"{len(to_predict)} to predict")


def store(args, settings):
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(args.fasta, "fasta")}
    results = read_nucimport_outputs(args.results)

    missing = [prot_id for prot_id in sequences if prot_id not in results]
    write_atomic(args.missing, [f"{prot_id}\n" for prot_id in missing])
    if missing:
        message = (f"{len(missing)} of {len(sequences)} proteins sent to NucImport got no result. "
                   f"They are listed in {args.missing}")
        if len(missing) > MAX_MISSING_FRACTION * len(sequences):
            logging.error(f"{message}. Please check the run_nucimport logs. Abort!")
            sys.exit(1)
        logging.warning(f"{message} and will be predicted again in the next run.")

    if args.no_cache:
        logging.info("NucImport cache turned off (nls_cache): no results stored")
        return

    # Another run may have stored some of these sequences in the meantime
    stored = read_store(args.cache_dir, settings)
    new_lines = []
    for prot_id, seq in sequences.items():
        seq_hash = sequence_hash(seq)
        if prot_id in results and seq_hash not in stored:
            new_lines.append("\t".join([seq_hash, *settings] + results[prot_id]) + "\n")
            stored[seq_hash] = results[prot_id]
    if not new_lines:
        logging.info("No new NucImport results to store")
        return

    os.makedirs(args.cache_dir, exist_ok=True)
    store_file = os.path.join(args.cache_dir, f"{args.name}_{time.strftime('%Y%m%d-%H%M%S')}_{os.getpid()}.tsv")
    write_atomic(store_file, ["\t".join(STORE_COLUMNS) + "\n"] + new_lines)
    logging.info(f"{len(new_lines)} new NucImport results stored in {store_file}")


def main():
    parser = argparse.ArgumentParser(description="Store of NucImport results, keyed by protein sequence")
    subparsers = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--fasta", required=True, help="Proteins of this run: the unique proteins for lookup, the proteins sent to NucImport for store")
    cache_group = common.add_mutually_exclusive_group(required=True)
    cache_group.add_argument("--cache_dir", help="Folder of the store files")
    cache_group.add_argument("--no_cache", action="store_true", help="Cache turned off: lookup sends every protein to NucImport and store only checks the results")
    common.add_argument("--model", required=True, help="NucImport training model, e.g. Model4")
    common.add_argument("--organism", required=True, help="NucImport organism option, Mouse or Yeast")
    common.add_argument("--id_flag", required=True, help="NucImport ID option, ID=F or ID=T")
    common.add_argument("--jar", required=True, help="NucImport jar, whose MD5 is part of the key")

    lookup_parser = subparsers.add_parser("lookup", parents=[common], help="Split the proteins into cached and to predict")
    lookup_parser.add_argument("--to_predict", required=True, help="Output FASTA with the proteins not in the cache")
    lookup_parser.add_argument("--cached", required=True, help="Output with the cached results, in NucImport format")

    store_parser = subparsers.add_parser("store", parents=[common], help="Add new NucImport results to the store")
    store_parser.add_argument("--results", required=True, nargs="+", help="NucImport output files")
    store_parser.add_argument("--name", required=True, help="Start of the store file name, e.g. Hsapiens_ensembl")
    store_parser.add_argument("--missing", required=True, help="Output list of the proteins without a NucImport result")

    args = parser.parse_args()
    # ID=T makes NucImport use the protein ID, so the result would no longer depend on the sequence only
    if args.id_flag.upper() != "ID=F" and not args.no_cache:
        logging.error(f"The NucImport cache needs ID=F, not {args.id_flag}. Abort!")
        sys.exit(1)

    settings = (args.model, args.organism, args.id_flag, file_md5(args.jar))
    if args.command == "lookup":
        lookup(args, settings)
    else:
        store(args, settings)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
