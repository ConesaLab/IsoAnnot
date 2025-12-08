#!/usr/bin/env python3
import argparse, sys, os, re, logging, traceback

def mirbase_parser(mirbase_file):
    """
    Parses the miRBase .dat file.
    Only saves products (mature miRNAs) that explicitly have the
    tag '/evidence=experimental'.
    """
    mirbase_dict = {}
    
    logging.info(f"Reading miRBase file: {mirbase_file}...")

    with open(mirbase_file, "r") as mirbase:
        line = mirbase.readline()
        
        name_sp = "Unknown"
        current_product = ""
        
        while line:
            id_line = line[0:2]
            content = line[2:].strip()

            # --- 1. DETECT SPECIES ---
            if id_line == "DE": 
                parts = content.split()
                if len(parts) >= 2:
                    name_sp = parts[0] + " " + parts[1]
                else:
                    name_sp = "Unknown"

            # --- 2. DETECT PRODUCTS AND EVIDENCE ---
            elif id_line == "FT": 
                if content.startswith("/pro"):
                    match = re.search(r'\"(.+)\"', content)
                    if match: 
                        current_product = match.group(1)
                
                elif content.startswith("/evi"):
                    # Only accept if it says "experimental" and NOT "not_experimental"
                    # This discards non-functional passenger strands (star strands)
                    if "experimental" in content and "not_experimental" not in content:
                        if name_sp not in mirbase_dict:
                            mirbase_dict[name_sp] = []
                        
                        if current_product != "":
                            mirbase_dict[name_sp].append(current_product)

            # --- 3. END OF ENTRY ---
            elif id_line == "//": 
                name_sp = "Unknown"
                current_product = ""

            line = mirbase.readline()
            
    return mirbase_dict

def mirwalkfilter(mirwalk_file, mirbase_dict, species, mirwalk_filtered, min_score=0.0):
    """
    Filters miRWalk based on:
    1. Existence of the miRNA in miRBase (Experimental).
    2. Score (binding probability) >= min_score.
    
    IMPROVEMENT: Always creates the output file to avoid errors in Snakemake.
    """
    
    # --- CRITICAL STEP: ALWAYS OPEN OUTPUT ---
    # This ensures the file exists (even if empty) so Snakemake does not fail.
    with open(mirwalk_filtered, "w") as mw_out:
        
        # 1. Check if species exists
        if species not in mirbase_dict:
            logging.warning(f"WARNING: Species '{species}' was NOT found in the miRBase file.")
            logging.warning("An EMPTY output file has been generated and the process will continue.")
            return # Exit successfully (without writing anything else)

        # 2. Load valid miRNAs for that species
        valid_mirnas = set()
        for mirna_name in mirbase_dict[species]:
            valid_mirnas.add(mirna_name)
            
        logging.info(f"Species '{species}' found. {len(valid_mirnas)} valid miRNAs loaded into memory.")

        # 3. Filter the input file
        logging.info(f"Filtering {mirwalk_file} (Score min: {min_score})...")
        
        count_in = 0
        count_kept = 0
        count_score_fail = 0
        
        with open(mirwalk_file, "r") as mw_in:
            for line in mw_in:
                # Keep headers
                if line.startswith("miRNA") or line.startswith("ID"): 
                    mw_out.write(line)
                    continue 
                
                cols = line.strip().split("\t")
                if len(cols) < 6: continue 
                
                # CURRENT MIRWALK FORMAT (according to your file):
                # Col 0: miRNA (hsa-let-7a-5p)
                # Col 5: Probability (0.84...)
                mirna_id = cols[0]
                
                try:
                    score = float(cols[5])
                except ValueError:
                    continue

                count_in += 1
                
                # --- FILTERING LOGIC ---
                if mirna_id in valid_mirnas:
                    if score >= min_score:
                        mw_out.write(line)
                        count_kept += 1
                    else:
                        count_score_fail += 1
                    
        logging.info(f"Process finished.")
        logging.info(f"Read: {count_in}")
        logging.info(f"Discarded due to low Score (<{min_score}): {count_score_fail}")
        logging.info(f"Kept: {count_kept}")

def main():
    parser = argparse.ArgumentParser(description="miRNA filter by experimental evidence and score (Robust)")
    
    parser.add_argument("--mirbase_file", required=True, help="miRBase .dat file")
    parser.add_argument("--mirwalk_file", required=True, help="miRWalk interactions file")
    parser.add_argument("--species", required=True, help="Scientific name of the species (e.g. 'Homo sapiens')")
    parser.add_argument("--mirwalk_output", required=True, help="Filtered output file")
    parser.add_argument("--score", type=float, default=0.0, help="Minimum binding probability (0.0 - 1.0). Default: 0")
    
    args = parser.parse_args()

    try:
        mirbase_dict = mirbase_parser(args.mirbase_file)
        mirwalkfilter(args.mirwalk_file, mirbase_dict, args.species, args.mirwalk_output, args.score)
        
    except Exception as ex:
        logging.error(f"Fatal error: {ex}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()