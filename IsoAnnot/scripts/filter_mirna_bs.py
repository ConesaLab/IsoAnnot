#!/usr/bin/env python3
import argparse, sys, os, re, logging, traceback

def mirbase_parser(mirbase_file, target_species):
    """
    Parses the miRBase .dat file.
    Only stores products for the 'target_species'.
    Ignores all other species.
    """
    valid_mirnas = set()
    
    logging.info(f"Reading miRBase file: {mirbase_file}...")
    logging.info(f"Extracting ONLY data for: '{target_species}'")

    with open(mirbase_file, "r") as mirbase:
        line = mirbase.readline()
        
        is_target_species = False
        current_product = ""
        
        while line:
            id_line = line[0:2]
            content = line[2:].strip()

            if id_line == "DE": 
                parts = content.split()
                if len(parts) >= 2:
                    name_sp = parts[0] + " " + parts[1]
                    
                    if name_sp == target_species:
                        is_target_species = True
                    else:
                        is_target_species = False
                else:
                    is_target_species = False

            elif id_line == "FT" and is_target_species: 
                if content.startswith("/pro"):
                    match = re.search(r'\"(.+)\"', content)
                    if match: 
                        current_product = match.group(1)
                
                elif content.startswith("/evi"):
                    if "experimental" in content and "not_experimental" not in content:
                        if current_product != "":
                            valid_mirnas.add(current_product)

            elif id_line == "//": 
                is_target_species = False
                current_product = ""

            line = mirbase.readline()
            
    return valid_mirnas

def mirwalkfilter(mirwalk_file, valid_mirnas, species, mirwalk_filtered, min_score=0.0):
    """
    Filters miRWalk.
    Input 'valid_mirnas' is a SET of IDs for the specific species.
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
                if line.startswith("miRNA") or line.startswith("ID") or line.startswith("#"): 
                    mw_out.write(line)
                    continue 
                
                cols = line.strip().split()
                if len(cols) < 6: continue 
                
                mirna_id = cols[0]
                
                try:
                    score = float(cols[5])
                except ValueError:
                    continue

                count_in += 1
                
                if mirna_id in valid_mirnas:
                    if score >= min_score:
                        mw_out.write(line)
                        count_kept += 1
                    else:
                        count_score_fail += 1
                    
        logging.info(f"Process finished.")
        logging.info(f"Binding sites read: {count_in}")
        logging.info(f"Binding sites discarded due to low Score (<{min_score}): {count_score_fail}")
        logging.info(f"Binding sites kept: {count_kept}")

def main():
    parser = argparse.ArgumentParser(description="miRNA filter (Optimized)")
    
    parser.add_argument("--mirbase_file", required=True, help="miRBase .dat file")
    parser.add_argument("--mirwalk_file", required=True, help="miRWalk interactions file")
    parser.add_argument("--species", required=True, help="Scientific name (e.g. 'Homo sapiens')")
    parser.add_argument("--mirwalk_output", required=True, help="Filtered output file")
    parser.add_argument("--score", type=float, default=0.0, help="Min binding probability")
    
    args = parser.parse_args()

    try:
        valid_mirnas_set = mirbase_parser(args.mirbase_file, args.species)
        
        mirwalkfilter(args.mirwalk_file, valid_mirnas_set, args.species, args.mirwalk_output, args.score)
        
    except Exception as ex:
        logging.error(f"Fatal error: {ex}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
