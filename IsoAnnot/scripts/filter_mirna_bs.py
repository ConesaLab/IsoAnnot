#!/usr/bin/env python3
"""
Author: Raúl López

Script to filter miRNA-target interactions from miRWalk using experimental evidence
from miRBase. It extracts valid miRNAs for a specific species and applies a
probability score threshold to the interaction data.
It also classifies each site by its seed match (TargetScan site types) and keeps only the
types asked for. The type is written as an extra last column.
"""

import argparse, sys, os, re, logging, traceback
from collections import Counter

WATSON_CRICK = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}
WOBBLE = {("G", "T"), ("T", "G")}

# Site types kept by each --site_types setting. 6mers have a marginal effect on average, and
# non-canonical sites (most of miRWalk's) bind without much repression (Agarwal et al. 2015, eLife).
ALL_SITE_TYPES = ["8mer", "7mer-m8", "7mer-A1", "6mer", "offset-6mer", "seed-GU", "noncanonical"]
SITE_TYPES = {
    "canonical": set(ALL_SITE_TYPES[:3]),
    "canonical_6mer": set(ALL_SITE_TYPES[:5]),
    "all": set(ALL_SITE_TYPES),
}

def duplex_partners(mirna_struct, target_struct):
    """
    Pairs the brackets of the miRWalk duplex structure.

    :return: for each miRNA position (0-based), the paired target position (0-based) or None.
    """
    structure = mirna_struct + target_struct
    n_mirna = len(mirna_struct)
    stack, pair = [], {}
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")":
            j = stack.pop()
            pair[i], pair[j] = j, i
    if stack:
        raise ValueError("unbalanced duplex structure")
    return [pair[i] - n_mirna if pair.get(i, -1) >= n_mirna else None for i in range(n_mirna)]

def seed_paired(mirna, target, partners, first, last, allowed):
    """True if miRNA positions first..last (1-based) pair with consecutive target bases."""
    start = partners[first - 1]
    if start is None:
        return False
    for k, pos in enumerate(range(first, last + 1)):
        j = partners[pos - 1]
        if j is None or j != start - k or (mirna[pos - 1], target[j]) not in allowed:
            return False
    return True

def site_type(duplex):
    """
    Seed match type of a miRWalk site, from its duplex "<miRNA>#<target>#<miRNA structure>#<target structure>"
    (both strands 5'->3'). 8mer: Watson-Crick pairs at miRNA positions 2-8 and an A in the target opposite
    position 1; 7mer-m8: 2-8; 7mer-A1: 2-7 and the A; 6mer: 2-7; offset-6mer: 3-8. seed-GU: 2-7 paired only
    when G:U wobbles count. Anything else is noncanonical.

    :return: the type, or None if the duplex cannot be read.
    """
    try:
        mirna, target, mirna_struct, target_struct = duplex.split("#")
        partners = duplex_partners(mirna_struct, target_struct)
    except (ValueError, IndexError):
        return None
    if len(mirna) != len(mirna_struct) or len(target) < len(target_struct) or len(mirna) < 8:
        return None

    # In about 20% of the sites the target structure is shorter than the target sequence: unpaired bases
    # are left out, nearly always at the 5' end. Use the offset where every pair is complementary.
    extra = len(target) - len(target_struct)
    for offset in [extra] + list(range(extra)):
        if all((mirna[i], target[j + offset]) in WATSON_CRICK | WOBBLE for i, j in enumerate(partners) if j is not None):
            break
    else:
        return None
    partners = [j + offset if j is not None else None for j in partners]

    # The A1 position is the target base right 3' of the one paired with miRNA position 2
    a1 = partners[1] is not None and partners[1] + 1 < len(target) and target[partners[1] + 1] == "A"
    if seed_paired(mirna, target, partners, 2, 8, WATSON_CRICK):
        return "8mer" if a1 else "7mer-m8"
    if seed_paired(mirna, target, partners, 2, 7, WATSON_CRICK):
        return "7mer-A1" if a1 else "6mer"
    if seed_paired(mirna, target, partners, 3, 8, WATSON_CRICK):
        return "offset-6mer"
    if seed_paired(mirna, target, partners, 2, 7, WATSON_CRICK | WOBBLE):
        return "seed-GU"
    return "noncanonical"

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

def mirwalkfilter(mirwalk_file, valid_mirnas, species, mirwalk_filtered, min_score=0.0, site_types="canonical"):
    """
    Filters miRWalk interaction files based on valid miRNA IDs, a score threshold and the seed match type.
    Kept lines get the site type as an extra last column.

    Args:
        mirwalk_file (str): Path to the raw miRWalk interactions file.
        valid_mirnas (set): Set of validated miRNA IDs from mirbase_parser.
        species (str): Target species name for logging purposes.
        mirwalk_filtered (str): Path to the output filtered file.
        min_score (float): Minimum binding probability (0.0 to 1.0).
        site_types (str): Key of SITE_TYPES: canonical, canonical_6mer or all.

    Returns:
        None
    """
    kept_types = SITE_TYPES[site_types]
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
        count_types = Counter()
        count_unreadable = 0

        with open(mirwalk_file, "r") as mw_in:
            for line in mw_in:
                # Preserve headers and comments
                if line.startswith("miRNA"):
                    mw_out.write(line.rstrip("\r\n") + "\tsite_type\n")
                    continue
                if line.startswith("ID") or line.startswith("#"):
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
                        kind = site_type(cols[4])
                        if kind is None:
                            count_unreadable += 1
                            kind = "noncanonical"
                        count_types[kind] += 1
                        if kind in kept_types:
                            mw_out.write(line.rstrip("\r\n") + f"\t{kind}\n")
                            count_kept += 1
                    else:
                        count_score_fail += 1

        logging.info(f"Process finished.")
        logging.info(f"Binding sites read: {count_in}")
        logging.info(f"Discarded due to low Score (<{min_score}): {count_score_fail}")
        for kind in ALL_SITE_TYPES:
            status = "kept" if kind in kept_types else "discarded"
            logging.info(f"Seed match {kind}: {count_types[kind]} ({status})")
        if count_unreadable:
            logging.warning(f"{count_unreadable} sites had a duplex that could not be read and were counted as noncanonical.")
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
    parser.add_argument("--site_types", choices=list(SITE_TYPES), default="canonical",
                        help="Seed match types to keep: canonical (8mer, 7mers), canonical_6mer (also 6mers) or all (default: %(default)s)")

    args = parser.parse_args()

    try:
        # Step 1: Parse miRBase to get the miRNAs for the species
        valid_mirnas_set = mirbase_parser(args.mirbase_file, args.species)

        # Step 2: Apply the filter to the miRWalk dataset
        mirwalkfilter(args.mirwalk_file, valid_mirnas_set, args.species, args.mirwalk_output, args.score, args.site_types)

    except Exception as ex:
        logging.error(f"Fatal error: {ex}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    # Configure logging to provide clean status updates to the terminal
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
