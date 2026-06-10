#!/usr/bin/env python3
"""
Author: Raúl López

Script to filter miRNA-target interactions from miRWalk using experimental evidence 
from miRBase. It extracts valid miRNAs for a specific species and applies a 
probability score threshold to the interaction data.
"""

import argparse, sys, os, re, logging, traceback

def mirbase_parser(mirbase_file, target_species):
    """
    Parses the miRBase .dat file to identify experimentally validated miRNAs.

    Args:
        mirbase_file (str): Path to the miRBase database file (.dat format).
        target_species (str): Scientific name of the species to filter (e.g., 'Homo sapiens').

    Returns:
        set: A set containing the names of miRNAs that have experimental evidence.
    """
    valid_mirnas = set()

    logging.info(f"Reading miRBase file: {mirbase_file}...")
    logging.info(f"Extracting ONLY data for: '{target_species}'")

    with open(mirbase_file, "r") as mirbase:
        line = mirbase.readline()

        is_target_species = False
        current_product = ""

        # miRBase .dat files are structured in blocks separated by '//'
        while line:
            id_line = line[0:2]
            content = line[2:].strip()

            # DE line contains the species description
            if id_line == "DE":
                parts = content.split()
                if len(parts) >= 2:
                    name_sp = parts[0] + " " + parts[1]
                    is_target_species = (name_sp == target_species)
                else:
                    is_target_species = False

            # FT lines contain features like products and evidence types
            elif id_line == "FT" and is_target_species:
                if content.startswith("/pro"):
                    # Extract product name inside quotes
                    match = re.search(r'\"(.+)\"', content)
                    if match:
                        current_product = match.group(1)

                elif content.startswith("/evi"):
                    # Only keep products with experimental evidence
                    if "experimental" in content and "not_experimental" not in content:
                        if current_product != "":
                            valid_mirnas.add(current_product)

            # Reset status at the end of each record block
            elif id_line == "//":
                is_target_species = False
                current_product = ""

            line = mirbase.readline()

    return valid_mirnas

def mirwalkfilter(mirwalk_file, valid_mirnas, species, mirwalk_filtered, min_score=0.0):
    """
    Filters miRWalk interaction files based on valid miRNA IDs and a score threshold.

    Args:
        mirwalk_file (str): Path to the raw miRWalk interactions file.
        valid_mirnas (set): Set of validated miRNA IDs from mirbase_parser.
        species (str): Target species name for logging purposes.
        mirwalk_filtered (str): Path to the output filtered file.
        min_score (float): Minimum binding probability (0.0 to 1.0).

    Returns:
        None
    """
    with open(mirwalk_filtered, "w") as mw_out:
        if not valid_mirnas:
            logging.warning(f"WARNING: No experimental miRNAs found for species '{species}' in miRBase.")
            logging.warning("An EMPTY output file has been generated.")
            return

        logging.info(f"Species '{species}' found. {len(valid_mirnas)} valid miRNAs loaded.")
        logging.info(f"Filtering {mirwalk_file} (Score min: {min_score})...")

        count_in = 0
        count_kept = 0
        count_score_fail = 0

        with open(mirwalk_file, "r") as mw_in:
            for line in mw_in:
                # Preserve headers and comments
                if line.startswith("miRNA") or line.startswith("ID") or line.startswith("#"):
                    mw_out.write(line)
                    continue

                cols = line.strip().split()
                if len(cols) < 6:
                    continue

                mirna_id = cols[0]
                try:
                    score = float(cols[5]) # Probability score is typically in the 6th column
                except ValueError:
                    continue

                count_in += 1

                # Filter by presence in validated set AND score threshold
                if mirna_id in valid_mirnas:
                    if score >= min_score:
                        mw_out.write(line)
                        count_kept += 1
                    else:
                        count_score_fail += 1

        logging.info(f"Process finished.")
        logging.info(f"Binding sites read: {count_in}")
        logging.info(f"Discarded due to low Score (<{min_score}): {count_score_fail}")
        logging.info(f"Binding sites kept: {count_kept}")

def main():
    """
    Command-line interface for the miRNA interaction filter.
    """
    parser = argparse.ArgumentParser(description="miRNA filter based on miRBase evidence and miRWalk scores")

    parser.add_argument("--mirbase_file", required=True, help="miRBase .dat file")
    parser.add_argument("--mirwalk_file", required=True, help="miRWalk interactions file")
    parser.add_argument("--species", required=True, help="Scientific name (e.g. 'Homo sapiens')")
    parser.add_argument("--mirwalk_output", required=True, help="Filtered output file")
    parser.add_argument("--score", type=float, default=0.0, help="Min binding probability (threshold)")

    args = parser.parse_args()

    try:
        # Step 1: Parse miRBase to get the miRNAs for the species
        valid_mirnas_set = mirbase_parser(args.mirbase_file, args.species)

        # Step 2: Apply the filter to the miRWalk dataset
        mirwalkfilter(args.mirwalk_file, valid_mirnas_set, args.species, args.mirwalk_output, args.score)

    except Exception as ex:
        logging.error(f"Fatal error: {ex}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    # Configure logging to provide clean status updates to the terminal
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
