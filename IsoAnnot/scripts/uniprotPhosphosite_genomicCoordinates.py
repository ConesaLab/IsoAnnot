#!/usr/bin/env python3
"""
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez
Modified by: Alessandra Martinez

Script that gets the genomic coordinates of uniprot proteins and their features.
"""

import argparse, sys, os, re, logging, csv, itertools, traceback
from operator import itemgetter
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from IsoAnnot import openfile, get_consecutive_parts, read_chr_ref_acc
from CoordinateMapper.CoordinateMapper import CoordinateMapper
from CoordinateMapper.MapPositions import ProteinPosition, ProteinPositionError


def get_fasta_sequences(fasta_files, matching_regex, keep_version="False", filter_function=None, reversed=False):
    '''
    Creates a dictionary that stores information from protein fasta files {ID:sequence}

    Args:
        Fasta files (list): uniprot, ensembl and refseq.
        Matching regex (str): regex used to extract sequence IDs.
        keep_version (boolean): keep version of sequence identifier in database
                              (ENSP3423143.1 -> ENSP3423143)
        filter function: function used to filter fasta records.
        reversed (bool): 
    
    Returns: 
        output results (dict): {ID:sequence}
    '''
    if reversed:
        output_results = defaultdict(list)
    else:
        output_results = {}

    # Regex used to capture sequence IDs
    capture_regex = re.compile(matching_regex)

    for fasta_file in fasta_files:
        logging.info(f"Reading fasta file {fasta_file}...\n")
        with openfile(fasta_file, "rt") as read_handler:
            for record in SimpleFastaParser(read_handler):
                # Filter record if function provided TODO??
                if filter_function and not filter_function(record[0]):
                    continue

                # Extract name from the FASTA identifier
                # Regex must contain only one group
                id_name = capture_regex.match(record[0]).group(1) 
                if keep_version.upper()=="FALSE":
                    id_name = id_name.split(".")[0]

                if not reversed:
                    output_results[id_name] = record[1]
                else:
                    output_results[record[1]].append(id_name)

    return(output_results)


def protein_match(seq1, seq2, max_mismatches=3):
    """
    Pairwise alignment of protein sequences allowing no 
    gaps and a maximum number of mismatches.

    Args:
        seq1 (str): protein sequence 1
        seq2 (str): protein sequence 2
        max_mismatches (int): maximum number of mismatches allowed

    Returns:
        Boolean indicating if sequences match
    """
    mismatch = 0

    if seq1 == None or seq2 == None:
        return False

    else: 
        if len(seq1) == len(seq2):
            for i in range(0, len(seq1)):
                if seq1[i] != seq2[i]:
                    mismatch +=1
            
            if mismatch <= max_mismatches:
                return True
        else:
            return False


def get_uniprot_associations(uniprot_file, proteins_associations, valid_features):
    '''
    Reads Uniprot Parsed file and gets the information regarding protein features an associations to other databases.

    Args:
        uniprot_file (str)
        proteins_associations (dict): dictionary {ID:sequence} with the proteins we want to inspect.
        valid_features (list): protein features we are interested in. 

    Returns:
        output_dict (dict): 
    '''
    output_dict = {
        "associations": defaultdict(lambda: defaultdict(lambda : {"transcripts": [], "proteins": []})),
        "features": defaultdict(list),
        "ids_all": set([]),
        "ids_with_database": set([]),
        "ids_with_match": set([]),
    }

    # Define the positions in the parsed uniprot file for the fields:
    # [referenceID, referenceTranscriptID]
    databases_indexes = {
        "RefSeq": [4, 5],
        "Ensembl": [5, 4],
        "EnsemblPlants": [5, 4],
        "EnsemblMetazoa": [5, 4]
    }

    logging.info(f"Analyzing uniprot data {uniprot_file}")

    with openfile(uniprot_file, "rt") as read_handler:
        tsv_reader = csv.reader(read_handler, delimiter="\t")
        # Row structure might be different depending on the row type and database 
        # But some fields are always the same
        for row in tsv_reader:
            uniprot_id, taxa_id_row, row_type, db_or_class, start_pos, end_pos = row[:6]  

            output_dict['ids_all'].add(uniprot_id)
            row_desc = row[6] if row[6:] else "NA"
            
            # Rows containing a database association (DR) with Ensembl/Refseq
            if row_type == "DR" and  db_or_class in databases_indexes:
                # Get protein and transcript IDs in those databases
                refID, refTransID = itemgetter(*databases_indexes.get(db_or_class))(row)
                output_dict['ids_with_database'].add(uniprot_id)
                uniprot_isoform = ""

                if db_or_class == "RefSeq":
                    if ". [" in refTransID:
                        uniprot_isoform = refTransID.split(". [")[1].strip("]")
                        refTransID = refTransID.split(". [")[0]
                    refID = refID.split(".")[0]
                    refTransID = refTransID.split(".")[0] 

                if db_or_class == "EnsemblPlants": #keep version
                    # refID = refID.split(".")[0]
                    # refTransID = refTransID.split(".")[0] 
                    refGeneID = row[-1] 
                    if ". [" in refGeneID:
                        uniprot_isoform = refGeneID.split(". [")[1].strip("]")

                if db_or_class == "Ensembl" or  db_or_class == "EnsemblMetazoa":
                    refID = refID.split(".")[0]
                    refTransID = refTransID.split(".")[0] 
                    refGeneID = row[-1] 
                    if ". [" in refGeneID:
                        uniprot_isoform = refGeneID.split(". [")[1].strip("]")

                if len(uniprot_isoform) == 0:
                    uniprot_isoform = uniprot_id

                # Keep the association if the UniprotKB protein matches the protein 
                # encoded by the associated Ensembl/RefSeq transcript
                # allowing a maximum of 3 mismatches 
                if refID in proteins_associations:
                    if protein_match(proteins_associations.get(uniprot_isoform), proteins_associations.get(refID)):
                        output_dict['ids_with_match'].add(uniprot_id)
                        output_dict['associations'][uniprot_id][db_or_class]['proteins'].append(refID)
                        output_dict['associations'][uniprot_id][db_or_class]['transcripts'].append(refTransID)
                    # canonical isoforms : have their version number in uniprot parsed file, but not in the fasta file
                    # eg: P26378-6
                    elif protein_match(proteins_associations.get(uniprot_id), proteins_associations.get(refID)):
                        output_dict['ids_with_match'].add(uniprot_id)
                        output_dict['associations'][uniprot_id][db_or_class]['proteins'].append(refID)
                        output_dict['associations'][uniprot_id][db_or_class]['transcripts'].append(refTransID)
            # For rows containing a protein FEATURE
            if row_type == "FT" and uniprot_id in output_dict["ids_with_match"] and \
                (valid_features is None or db_or_class in valid_features) and start_pos.isdigit() and end_pos.isdigit():
                feature_instance = SeqFeature(FeatureLocation(int(start_pos), int(end_pos)),
                                                type=db_or_class,
                                                id=uniprot_id,
                                                qualifiers={"desc": row_desc})
                output_dict['features'][uniprot_id].append(feature_instance)
    
    return(output_dict)


def get_cds_from_gtf(gtf_file, chrom_ref={}, biomart_host = ""):
    '''
    Gets CDS coordinates from GTF files 

    Args:
        gtf_file (str): Ensembl or Refseq.
        match_file (str): chromosome accession file needed for Refseq database.

    Returns:
        feature_cds (dict): 
    '''
    logging.info(f"Reading GTF file {gtf_file}")

    protein_cds = defaultdict(list)
    feature_cds = {}
    strand_conversion = {"+": 1, "-": -1}
    
    with open(gtf_file, "r") as gtf_content:
        for line in gtf_content:
            if line[0] != "#":
                if line.split("\t")[2]=="CDS":
                    try:
                        protein_id = line.split("protein_id")[1].split(";")[0].replace('"', '').strip()
                    except IndexError:
                        continue
                    if "." in protein_id and not len(chrom_ref) == 0 and "plants" not in biomart_host:
                        protein_version = protein_id.split(".")[1]
                        protein_id = protein_id.split(".")[0]

                    chr = line.split("\t")[0]
                    chr_name = chrom_ref.get(chr, chr)
                    strand = line.split("\t")[6]
                    biopython_strand = strand_conversion.get(strand, strand)
                    start=int(line.split("\t")[3])
                    end=int(line.split("\t")[4])
                    seq_record = SeqFeature(FeatureLocation(start -1, end,
                                                biopython_strand), id=chr_name)
                    protein_cds[protein_id].append(seq_record)

    for protein_id, feature_list in protein_cds.items():
        location_list =  [feature.location for feature in feature_list] 

        # CompoundLocation does not accept less than 2 parts
        if len(location_list) > 1:
            location_object = CompoundLocation(location_list)
        else:
            location_object = location_list.pop()

        feature_cds[protein_id] = SeqFeature(location_object, type="CDS",
                                             id=feature_list[0].id, qualifiers={"ref": protein_id})
                                            # id = chr

    return(feature_cds)


def get_phosphosite_modifications(filename):
    '''
    Parses Phosphosite files and creates a dictionary with information about protein modifications

    Args:
        filename (str)

    Returns:
        modifications_dict (dict)
    '''
    modifications_dict = defaultdict(list)
    capture_regex = re.compile("(\D+)(\d+)\-(.+)")

    mod_to_feature = {
        "ga": "CARBOHYD",
        "gl": "CARBOHYD"
    }

    mod_to_desc = {
        "ga": "O-linked (GalNAc...)",
        "gl": "O-linked (GlcNAc)",
        "sm": "Sumoylation",
        "ub": "Ubiquitination",
        "me": "methyl",
        "m1": "methyl",
        "m2": "dimethyl",
        "m3": "trimethyl",
        "p": "phospho",
        "ac": "acetyl"
    }

    aa_to_name = {"C": "cysteine", "D": "aspartate", "S": "serine", "Q": "glutamine", "K": "lysine",
                  "I": "isoleucine", "P": "proline", "T": "threonine", "F": "phenylalanine", "N": "asparagine",
                  "G": "glycine", "H": "histidine", "L": "leucine", "R": "arginine", "W": "tryptophan",
                  "A": "alanine", "V": "valine", "E": "glutamate", "Y": "tyrosine", "M": "methionine"}

    with openfile(filename, "rt") as read_handler:
        tsv_reader = csv.reader(itertools.islice(read_handler, 4, None), delimiter="\t")
        # Get modification information
        for phosphosite_row in tsv_reader:
            uniprot_id = phosphosite_row[2]
            phospho_mod = phosphosite_row[4]
            (motif, aa_position, mod_type) = capture_regex.match(phospho_mod).group(1, 2, 3)
            desc = mod_to_desc.get(mod_type, "")

            if mod_type in ["p", "m1", "m2", "m3", "m", "ac"]:
                desc += aa_to_name.get(motif)

            # Dictionary {proteinID:[list_of_modifications]}
            modifications_dict[uniprot_id].append({"motif": motif,
                                                   "modification": mod_type,
                                                   "aa_position": aa_position,
                                                   "feature": mod_to_feature.get(mod_type, "MOD_RES"),
                                                   "desc": desc,
                                                   "mod": phospho_mod})

    return(modifications_dict)


def main():
    try:
        parser = argparse.ArgumentParser(description="UniprotPhosphosite GenomicCoordinates")
        parser.add_argument("--uniprot_fasta", nargs="+", required=True)
        parser.add_argument("--refseq_fasta", nargs="+", required=False, default=[])
        parser.add_argument("--ensembl_fasta", nargs=1, required=True)
        parser.add_argument("--output_protein", required=True)
        parser.add_argument("--output_domain", required=True)
        parser.add_argument("--phosphosite_files", nargs="+", required=True)
        parser.add_argument("--uniprot_parsed", required=True)
        parser.add_argument("--uniprot_valid_features", 
                            default=["ACT_SITE", "BINDING", "CARBOHYD", "COILED", "COMPBIAS", "CROSSLNK",
                                     "DISULFID", "DNA_BIND", "INTRAMEM", "LIPID", "METAL", "MOD_RES", "MOTIF",
                                     "REGION", "SITE", "TRANSMEM", "ZN_FING"])
                           
        parser.add_argument("--refseq_gtf", required=False, default=None)
        parser.add_argument("--ensembl_gtf", required=True)
        parser.add_argument("--chr_ref", required=False, default=None)
        parser.add_argument('--biomart_host', nargs="?", const="http://www.ensembl.org")

        args = parser.parse_args()

        if "plants" not in args.biomart_host:
            keep_version_ensembl = "FALSE"
            keep_version_refseq = "FALSE"
        else:
            keep_version_ensembl = "TRUE" #In EnsemblPlants version differenciates different transcripts
            keep_version_refseq = "FALSE"

        # Refseq chromosomes
        if args.chr_ref:
            logging.info(f"Reading chr accession table {args.chr_ref}")
            refseqChrom = read_chr_ref_acc(args.chr_ref)
            logging.info(f"chr_conversion {refseqChrom}")

        # Compile 3 databases in one dictionary containing protein identifiers and the sequence.
        uniprot_proteins = get_fasta_sequences(fasta_files=args.uniprot_fasta,
                                                 matching_regex=".*?\\|(.*?)\\|",
                                                 keep_version="FALSE")
        logging.info(f"Len of uniprot proteins: {len(uniprot_proteins)}")

        refseq_proteins = get_fasta_sequences(fasta_files=args.refseq_fasta,
                                                matching_regex="(.*?)\\.?\\s",
                                                keep_version=keep_version_refseq)
        logging.info(f"Len of refseq proteins: {len(refseq_proteins)}")

        ensembl_proteins = get_fasta_sequences(fasta_files=args.ensembl_fasta,
                                                 matching_regex="(.*?)\\.?\\s",
                                                 keep_version=keep_version_ensembl) 
        logging.info(f"Len of ensembl proteins: {len(ensembl_proteins)}")

        proteins_dict = {**uniprot_proteins, **refseq_proteins, **ensembl_proteins}

        # Get links to other databases.
        uniprot_associations = get_uniprot_associations(uniprot_file=args.uniprot_parsed,
                                                        proteins_associations=proteins_dict,
                                                        valid_features=args.uniprot_valid_features)
        logging.info(f"Len of uniprot associations: {len(uniprot_associations['associations'])}")
        logging.info(f"Len of uniprot features: {len(uniprot_associations['features'])}")

        # TODO: write a "match" and "noMatch" files if required

        # Get CDS information from databases
        # TODO: make this handler more general
        cds_ensembl = get_cds_from_gtf(args.ensembl_gtf, biomart_host=args.biomart_host)

        cds_refseq={}
        if args.refseq_gtf:
            cds_refseq = get_cds_from_gtf(args.refseq_gtf,
                                        chrom_ref=refseqChrom)  
    
        logging.info(f"Len of ensembl CDS: {len(cds_ensembl)}")
        logging.info(f"Len of refseq CDS: {len(cds_refseq)}")
        
        gtf_cds = {**cds_ensembl, **cds_refseq}

        uniprot_cds = defaultdict(list)

        strand_back_conversion = {1: "+", -1: "-"}
	
        def _retrieve_parts_positions(feature_gstart, feature_gend, exon_list):
            '''
            For a specific feature of a protein, it gets all the genomic 
            positions that, when translated, give rise to aminoacids are part 
            of the feature.

            Args:
                Protein_start (int): genomic start position of protein feature
                Protein_end (int): genomic end position of protein feature
                Exon_list (list): list of genomic positions that are part of 
                                  CDS.
            Returns:
                output_dict (dict): dictionary with genomic start and end
                genomic positions of the feature and its length.
            '''
            # Back to 1-based coordinates
            feature_gstart += 1 
            exon_list = [i+1 for i in exon_list]

            # Get genomic positions of the protein feature
            domain_gpos = sorted([feature_gstart, feature_gend])
            domain_gpos = set(range(domain_gpos[0], domain_gpos[1] + 1))

            # Intersection with mapper exon list
            common_pos = domain_gpos & set(exon_list)

            output_dict = {
                "start": [],
                "end": [],
                "len": 0
            }

            for loc_range in get_consecutive_parts(common_pos):
                output_dict["start"].append(loc_range[0])
                output_dict["end"].append(loc_range[-1])
                output_dict["len"] += loc_range[-1] - loc_range[0] + 1

            return output_dict


        with open(args.output_protein, "w") as output_protein:

            # OUTPUT 1: Protein file stores protein level information
            output_protein_tsv = csv.writer(output_protein, delimiter="\t")

            # Search associated protein and CDS in GTF files
            for uniprot_id, item_dbs in uniprot_associations['associations'].items():
                # Search in GTFs
                for db_id, db_data in item_dbs.items():
                    for protein_id in db_data.get("proteins"):
                        # Default values when no match is found
                        protein_len = len(proteins_dict.get(protein_id, 0))
                        protein_chr = protein_strand = cds_len = "NA"
                        match_tag = "NOT_id"
                        
                        # If protein CDS is found but does not match
                        if protein_id in gtf_cds:
                            protein_cds = gtf_cds.get(protein_id) 

                            # Override values
                            match_tag = "CDS_incomplete"
                            cds_len = len(protein_cds)
                            protein_chr = protein_cds.id
                            protein_strand = protein_cds.location.strand
                            protein_strand = strand_back_conversion.get(protein_strand, None)

                            # Complete the SeqFeature attributes if CDS length corresponds to protein length
                            if cds_len == (protein_len * 3):  
                                # Complete the SeqFeature attributes
                                # protein_cds.seq = proteins_dict.get(uniprot_id) #REVSIAR
                                protein_cds.qualifiers["uniprot"] = uniprot_id
                                #Store all the transcripts (CDSs) associated to a uniprot entry (main and secondary isoforms)
                                uniprot_cds[uniprot_id].append(protein_cds)
                                match_tag = "ok"

                            output_protein_tsv.writerow([
                                uniprot_id,     # Protein Uniprot ID
                                protein_id,     # Associated Protein ID in other database
                                protein_chr,    # Chromosome
                                protein_strand, # Strand
                                protein_len,    # Protein length (aa)
                                cds_len,        # Protein lenght (nt)
                                # Back to 1-based coordinates
                                list(map(lambda x: int(x.start) +1, protein_cds.location.parts)), # Start position of each CDS part
                                list(map(lambda x: int(x.end), protein_cds.location.parts)),  # End position of each CDS part
                                match_tag
                            ])

            # OUTPUT 2: Domain file contains info about domains, features and modifications in each protein 
            with  open(args.output_domain, "w") as output_domain:
                output_domain_tsv = csv.writer(output_domain, delimiter="\t")
                
                for uniprot_id, features in uniprot_associations['features'].items():
                    protein_sequence = proteins_dict.get(uniprot_id, None)  

                    # Get feature position in protein   
                    for prot_feature in features:
                        # Change 0-based half-open to 0-based closed coordinates needed by CoordinateMapper
                        feature_pstart = ProteinPosition(prot_feature.location.start, index=None)
                        feature_pend = ProteinPosition(prot_feature.location.end - 1, index=None)

                        # Get genomic coordinates of feature for all the isoforms associated to a uniprot protein
                        for cds in uniprot_cds.get(uniprot_id, []):
                            cm = CoordinateMapper(cds)

                        # Get genomic coordinates of the feature (0-based half open)
                            try:
                                feature_gstart = int(cm.p2g(feature_pstart)[0])
                                feature_gend = int(cm.p2g(feature_pend)[1]) +1

                            except ProteinPositionError:
                                continue
                            
                            parts_in_cds =  _retrieve_parts_positions(feature_gstart, feature_gend, cm.exon_list)
                            motif_sequence = protein_sequence[prot_feature.location.start:prot_feature.location.end]

                            output_domain_tsv.writerow([
                                prot_feature.type, # Feature type
                                prot_feature.qualifiers['desc'], # Description 
                                cds.id, # Chromosome
                                strand_back_conversion.get(cds.location.strand, None), # Strand
                                parts_in_cds["start"][0],   # Feature genomic start
                                parts_in_cds["end"][-1],    # Feature genomic end
                                ",".join(map(str, parts_in_cds["start"])), # Start position of feature in each CDS part
                                ",".join(map(str, parts_in_cds["end"])), # End position offeature in each CDS part
                                abs(feature_pend.pos - feature_pstart.pos) + 1, # Feature length (aa)
                                parts_in_cds["len"],    # Feature length (nt)
                                motif_sequence,         # Feature sequence
                                "SwissProt",            # Source
                                uniprot_id,             # Protein Uniprot ID
                                cds.qualifiers['ref']   # Associated Protein ID in other database
                            ])

                # Phosphosite modifications
                for filename in args.phosphosite_files:
                    for uniprot_id, phospho_mods in get_phosphosite_modifications(filename).items():
                        protein_sequence = proteins_dict.get(uniprot_id, None)
                        
                        if protein_sequence:
                            # Get modification position in protein in 0-based coordinates
                            for phospho_mod in phospho_mods:
                                mod_pposition = int(phospho_mod.get("aa_position")) -1
                                # Check if sequence in protein matches to sequence in Phosphosite file                
                                if mod_pposition <= len(protein_sequence) -1:
                                    motif_pseq = protein_sequence[mod_pposition]

                                    for cds in uniprot_cds.get(uniprot_id, []):
                                        if motif_pseq == phospho_mod.get("motif"):
                                            cm = CoordinateMapper(cds)

                                            # Get genomic coordinates 0-based half-open
                                            try: 
                                                domain_ppos = cm.p2g(mod_pposition) 
                                            except ProteinPositionError:
                                                continue
                                            
                                            domain_ppos_start = int(domain_ppos[0])
                                            domain_ppos_end = int(domain_ppos[1]) +1
                                            parts_in_cds = _retrieve_parts_positions(domain_ppos_start, domain_ppos_end, cm.exon_list)
                                            
                                            output_domain_tsv.writerow([
                                            phospho_mod['feature'], # Feature type
                                            phospho_mod['desc'],    # Description
                                            cds.id,  # Chromosome
                                            strand_back_conversion.get(cds.location.strand, None), # Strand
                                            # 1-based coordinates
                                            domain_ppos_start +1,  # Feature genomic start
                                            domain_ppos_end,    # Feature genomic end
                                            ",".join(map(str, parts_in_cds["start"])),   # Start position of feature in each CDS part
                                            ",".join(map(str, parts_in_cds["end"])),     # End position offeature in each CDS part
                                            1,      # Feature length (aa)
                                            parts_in_cds["len"],   # Feature length (nt)
                                            motif_pseq,         # Feature sequence
                                            "PhosphositePlus",     # Source
                                            uniprot_id,            # Protein Uniprot ID
                                            cds.qualifiers['ref'], # Associated Protein ID in other database
                                            phospho_mod['mod']      # Modification
                                    ])  

                                                
    except Exception as ex:
        logging.error(str(ex))
        traceback.print_exc()
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
