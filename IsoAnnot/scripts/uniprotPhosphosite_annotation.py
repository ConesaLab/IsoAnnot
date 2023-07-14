#!/usr/bin/env python3
"""
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez 
Modified by: Alessandra Martinez

Script that projects protein features genomic coordinates into query isoforms and returns their protein coordinates.
""" 

#TODO avoid repeating identical annotations when the only difference is 
# the source Uniprot or Phosphosite

import argparse, sys, os, re, logging, csv, traceback
import pandas as pd
from collections import defaultdict
from IsoAnnot import openfile, strand_table, get_consecutive_parts, argparse_nullable
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from uniprotPhosphosite_genomicCoordinates import get_fasta_sequences
from CoordinateMapper.CoordinateMapper import CoordinateMapper


def get_coding_transcripts(classification_filename, protein_sequences):
    """
    Gets transcripts location for transcripts for which we have a protein sequence.

    Args:
        -classification_filename (str)
        -protein_sequences (dict)

    Returns:
        -transcript_dict (dict)

    """
    logging.info(f"Reading SQANTI classification file {classification_filename}")

    transcript_dict = {}

    with openfile(classification_filename, "rt") as read_handler:
        tsv_reader = csv.DictReader(read_handler, delimiter="\t")
        for line in tsv_reader:
            if line.get('coding') == "coding":
                transcript_id = line.get('isoform') # or associatedTranscript?
                # In mytranscripts there is no refProt, use the transcript id
                protein_id = line.get('refProt', transcript_id)

                # SQANTI classification file is 1-based, change to 0-based
                cds_pos = [int(line.get('CDS_start'))-1, (int(line.get('CDS_end')))]
                protein_seq = protein_sequences.get(protein_id, None)

                # Create  SeqFeature object for the transcript with CDS position, transcript ID,
                # protein ID and protein sequence.
                if protein_seq:
                    seq_transcript = SeqFeature(FeatureLocation(min(cds_pos), max(cds_pos)),
                                                type="transcript")
                    record_transcript = SeqRecord(Seq(protein_seq),
                                id=transcript_id,
                                name=protein_id,
                                features=[seq_transcript])
                    
                    # TODO: should this be a list instead?
                    transcript_dict[transcript_id] = record_transcript

    return(transcript_dict)


def get_genepred_exons(genepred_filename, match_file = None):
    """
    Get exons location info for each transcript from genepred file
    Args: 
        -genepred_filename (str)
        - chr_assoc (dict): Refseq chromosome accession
    
    Returns:
        -output_dict (dict)
    """
    chr_conversion = {}

    # Get chr from Refseq
    if match_file is not None:
        logging.info(f"Reading chr accession table {match_file}")
        with openfile(match_file, "rt") as read_handler:
            tsv_reader = csv.reader(filter(lambda row: row[0] != "#", read_handler), delimiter="\t")
            chr_conversion = {row[1]: row[0] for row in tsv_reader}


    logging.info(f"Reading GenePred file {genepred_filename}")
    output_dict = {}

    with openfile(genepred_filename, "rt") as read_handler:
        tsv_reader =  csv.reader(filter(lambda row: row[0] != "#", read_handler), delimiter="\t")

        for line in tsv_reader:
            transcript_id, chrom, strand = line[0:3]

            chrom = chr_conversion.get(chrom, chrom) 

            exons_start = line[8].split(",")[:-1] 
            exons_end = line[9].split(",")[:-1]

            # Exclude transcript without strand information
            if len(strand):
                strand = strand_table.get(strand, None)
                # Create a FeatureLocation object for each exon of the transcript
                # Change to 1-based coordinates and end+1 to include the last position 
                transcript_locations = [FeatureLocation(int(e_start)+1, int(e_end)+1, strand) for e_start, e_end in zip(exons_start, exons_end)]
                location_object = sum(transcript_locations)
                output_dict[transcript_id] = {"chr": chrom, "location": location_object}

    return(output_dict)


def get_uniprot_domains(domaingenomic_filename):
    """
    Get uniprot domain/feature info (1-based file)
    
    Args:
        -domaingenomic_filename (str)
    
    Output:
        -sorted_dict (dict)

    """
    logging.info(f"Reading uniprot data file {domaingenomic_filename}")
    output_dict = defaultdict(list)

    motif_regex = re.compile("^([^;\.{]+)") 

    allowed_features = ["ACT_SITE", "BINDING", "CA_BIND", "DNA_BIND", "MOD_RES", "CARBOHYD", "LIPID",
                        "DISULFID", "METAL", "NONT_STD", "COILED", "INTRAMEM", "MOTIF", "REGION", "SITE",
                        "COMPBIAS", "CROSSLNK", "TRANSMEM", "ZN_FING"]

    motif_feature_assign = {
        "DNA_BIND": "BINDING",
        "ZN_FING": "BINDING",
        "MOD_RES": "PTM",
        "CARBOHYD": "PTM",
        "LIPID": "PTM",
        "DISULFID": "PTM",
        "METAL": "BINDING",
        "CROSSLNK": "BINDING",
        "REGION": "MOTIF",
        "SITE": "MOTIF"
    }


    def _check_mod_res(line_contents):
        """
        Determine Post Transcriptional Modification type
        """
        search_col = line_contents[1].lower()
        output = "PTM_other"

        if any(mod in search_col for mod in ["phosphoserine", "phosphothreonine", "phosphotyrosine"]):
            output = motif_regex.match(search_col).group(1)
        elif "phospho" in search_col:
            output = "Phosphorilation_nonT_nonS_nonY"
        elif "acetyl" in search_col:
            output = "Acetylation"
        elif "methyl" in search_col:
            output = "Methylation"

        return(output)


    def _check_mod_comp(line_contents):
        """
        Determine comp bias residue
        
        """
        search_col = line_contents[1].lower()
        output = "Comp_bias_Other"

        if any(mod in search_col for mod in ["ser-rich", "Poly-ser"]):
            output = "Ser-comp_bias"
        elif any(mod in search_col for mod in ["pro-rich", "Poly-pro"]):
            output = "Pro-comp_bias"
        elif any(mod in search_col for mod in ["glu-rich", "Poly-glu"]):
            output = "Glu-comp_bias"

        return(output)


    def _check_mod_motif(line_contents):
        """
        Determine motif type
        
        """
        search_col = line_contents[1].lower()
        output = "Motif"

        if "nuclear_export" in search_col:
            output = "NES"
        elif any(mod in search_col for mod in ["nuclear localization", "nuclear import"]):
            output = "NLS"

        return(output)


    motif_id_assign = {
        "MOD_RES": _check_mod_res,
        "COMPBIAS": _check_mod_comp,
        "MOTIF": _check_mod_motif
    }

    motif_name_assign = {
        "COILED": "Coiled_coil",
        "INTRAMEM": "Intramembrane region",
        "NLS": "Nuclear Localization Signal",
        "NES": "Nuclear Export Signal"
    }

    motif_desc_assign = {
        "COILED": "Regions of coiled coil within the protein"
    }

    # Get motif info
    with openfile(domaingenomic_filename, "rt") as read_handler:
        tsv_reader =  csv.reader(filter(lambda row: row[0] != "#", read_handler), delimiter="\t")

        for line in (line for line in tsv_reader if line[0] in allowed_features):
            motif_type = line[0]
            database = line[11] 
            protein_uniprot = line[12]
            protein_id = line[13]

            # Motif ID is a requirement of other parameters, resolve it first.
            motif_id = motif_id_assign.get(motif_type, motif_type)

            if callable(motif_id):
                motif_id = motif_id(line)

            motif_annotations = {
                'id': motif_id.title(),
                'feature': motif_feature_assign.get(motif_type, motif_type),
                'desc': None,
                'name': motif_name_assign.get(motif_id, motif_id).title(),
                'database': database
            }

            try:
                motif_annotations['desc'] = motif_desc_assign.get(motif_id, motif_regex.match(line[1]).group(1).title())
            except:
                motif_annotations['desc'] = line[1]

            # Call the function if required
            motif_annotations = {key: (value(line) if callable(value) else value) for key, value in motif_annotations.items()}

            # Construct location object
            chrom_id = line[2]

            # TODO: keep compatibility with old format, multiple start/end positions, or remove it?
            # motif_start = int(int_extract_re.match(line[6]).group(1))
            # motif_end = int(int_extract_re.match(line[7]).group(1))

            motif_info = {
                "chrom_id": line[2],
                "seq_strand": strand_table.get(line[3], None), 
                "aa_sequence": line[10],
                "start_pos": line[6].split(","),
                "end_pos": line[7].split(","),
                "annotations": motif_annotations 
            }

            # Motif file is 1-based, keep it and add end +1 to include the last position
            try:

                motif_locs = [FeatureLocation(int(start_pos), int(end_pos) +1, \
                    motif_info.get("seq_strand")) for start_pos, end_pos in \
                    zip(motif_info.get("start_pos"), motif_info.get("end_pos"))]

                motif_seqfeature = SeqFeature(sum(motif_locs), type="motif")
                motif_record = SeqRecord(Seq(motif_info.get("aa_sequence")),
                                    id=protein_uniprot,
                                    name=protein_id,
                                    features=[motif_seqfeature],
                                    annotations=motif_info.get("annotations"))
            except ValueError as e:
                logging.error(e)
                continue

            # Group motifs by chromosome
            # TODO check if it's better to use protein_id as KEY to make the search faster
            output_dict[chrom_id].append(motif_record)

    sorted_dict = {k: sorted(v, key=lambda x: x.features[0].location.start) for k, v in output_dict.items()}

    return(sorted_dict)

# TODO check if this works for EnsemblPlants mytranscripts
def read_protein_association(filename):
    protein_assoc_dict = defaultdict(lambda : {"uniprot": [], "ensembl": [], "refseq": []})

    with open(filename, "r") as infile:
        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            if row[3] == "sequence":
                if row[2] == "trembl" or row[2] == "swissprot":
                    protein_assoc_dict[row[0]]["uniprot"].append(row[1])
                elif row[2] == "ensembl":
                    protein_assoc_dict[row[0]]["ensembl"].append(row[1])
                elif row[2] == "refseq":
                    if "." not in row[1]:                      
                        protein_assoc_dict[row[0]]["refseq"].append(row[1])
                    else: 
                        protein_assoc_dict[row[0]]["refseq"].append(row[1].split(".")[0])
    return protein_assoc_dict


def main():
    try:
        parser = argparse.ArgumentParser(description='UniprotPhosphosite')
        parser.add_argument('--orf_fasta', nargs=1, required=True) 
        parser.add_argument('--classification_file', nargs=None, required=True)
        parser.add_argument('--genepred_file', nargs=None, required=True)
        parser.add_argument('--uniprotmotif_file', nargs=None, required=True)
        parser.add_argument("--chr_ref", required=False, nargs='?', const=None, type=argparse_nullable)
        parser.add_argument('--output', nargs=None, required=True)
        parser.add_argument('--keep_version', nargs=None, default="False")
        parser.add_argument('--db', nargs=None, required=True)
        parser.add_argument('--protein_association', nargs='?', default=False, required=False, type=argparse_nullable)
        parser.add_argument('--biomart_host', nargs="?", const="http://www.ensembl.org")


        args = parser.parse_args()

        # Get info from files 
        protein_sequences = get_fasta_sequences(fasta_files=args.orf_fasta,
                                                matching_regex="(.*?)\\.?\\s",
                                                keep_version=args.keep_version)
                
        coding_transcripts = get_coding_transcripts(classification_filename=args.classification_file,
                                                    protein_sequences=protein_sequences)

        transcript_pred_exons = get_genepred_exons(genepred_filename=args.genepred_file, match_file=args.chr_ref)
        chrom_domains = get_uniprot_domains(domaingenomic_filename=args.uniprotmotif_file)

        protein_assoc = None
        if args.protein_association:
            protein_assoc = read_protein_association(args.protein_association)

        # Start 
        with open(args.output, "w") as output_file:
            logging.info(f"Starting coordinate transference process in {args.output}")
            output_file_tsv = csv.writer(output_file, delimiter="\t")

            # For each transcript, get exons and CDS genomic location
            for transcript_id, transcript_record in coding_transcripts.items():
                transcript_exons = transcript_pred_exons.get(transcript_id, None)
                # Avoid entering the loop if there is no info about the CHR  
                if transcript_exons and transcript_exons.get('chr', None) in chrom_domains:
                    # Get CDS location related to the transcript positions
                    transcript_location = transcript_record.features[0].location 
                    transcript_exons_loc = transcript_exons.get('location')
                    # Transform the transcript coordinates to CDS to use CoordinateMapper
                    cds_start_offset, cds_end_offset = ((int(transcript_location.start)), (int(transcript_location.end)))
                    all_transcript_positions = sorted(list(transcript_exons_loc), reverse=transcript_exons_loc.strand < 0)
                    # Select the CDS section (stop codon not included)
                    cds_positions = sorted(all_transcript_positions[(cds_start_offset):(cds_end_offset)])
                    # Create the new location objects by grouping them in continuous ranges
                    cds_locations = []
                    # Feature Locations modified to correctly create the coordinate mapper
                    # when sorting cds_loc_position, it  keeps the first position but removes the last one  
                    # therefore cds_locations have to be 1-based start, and end+1
                    for loc_range in get_consecutive_parts(cds_positions):
                        cds_locations.append(FeatureLocation(loc_range[0], loc_range[-1] +1,  
                                                             strand=transcript_exons_loc.strand))

                    cds_loc_position = sum(cds_locations)
                    cds_cm = CoordinateMapper(SeqFeature(cds_loc_position, type="CDS"))              

                    # For each valid transcript, look for motifs that are in the same chromosome as the transcript
                    # And get motif information if it corresponds to the same protein_id associated to the transcript
                    for motif_record in chrom_domains.get(transcript_exons.get('chr'), {}):
                        motif = motif_record.features[0]
                        motif_set = set(motif.location)


                        if args.db.upper() == "MYTRANSCRIPTS":
                            # For those transcripts whose products match known uniprot proteins
                            if transcript_id in protein_assoc and len(protein_assoc[transcript_id]["uniprot"]) > 0:
                                motif_protein_id = motif_record.name
                                if motif_protein_id in protein_assoc[transcript_id]["ensembl"] or motif_protein_id in protein_assoc[transcript_id]["refseq"]:
                                    if motif.location.start > cds_loc_position.end:
                                            break

                                    # TODO check if conditional could be restricted to issubset() only, since it should include all the other conditions
                                    # if motif.strand == transcript_exons_loc.strand and cds_positions[0] < motif.location.end and cds_positions[-1] > motif.location.start and motif_set.issubset(cds_positions):
                                    if motif.strand == transcript_exons_loc.strand and transcript_exons_loc.start < motif.location.end and transcript_exons_loc.end > motif.location.start and motif_set.issubset(set(cds_loc_position)): # linea antigua
                                        # Until now all motif positions could be included inside the CDS, but the later could
                                        # contain a in-between fragment between them, so we perform a reverse search.
                                        filtered_cds = set(filter(lambda x: motif.location.start <= x < motif.location.end, sorted(cds_loc_position)))
                                        # Check that the motif is inside the transcript and has only one part (consecutive sequences)
                                        if filtered_cds.issubset(motif_set):
                                            try:    
                                                if motif.location.strand == 1:
                                                    first_aa_position = int(cds_cm.g2p(motif.location.start))
                                                else:
                                                # For reverse strand,  take into account stop codon
                                                    first_aa_position = int(cds_cm.g2p(motif.location.end-3))
                            
                                                last_aa_position = first_aa_position + len(motif_record.seq)
                                            except:
                                                continue

                                    # Extract the motif sequence and, if it matches the one in the input file, annotate it
                                    transcript_sequence = transcript_record.seq[first_aa_position:last_aa_position]
                                    if motif_record.seq == transcript_sequence:
                                        annotations = motif_record.annotations
                                        gtf_attributes = f'ID={annotations["id"]}; Name={annotations["name"]}; Desc={annotations["desc"]}_{annotations["database"]}'

                                        output_file_tsv.writerow([
                                            transcript_id, #seqname
                                            "UniProtKB/Swiss-Prot_Phosphosite", # source
                                            motif_record.annotations.get('feature'), # feature
                                            first_aa_position + 1, # start (1-based)
                                            last_aa_position, # end
                                            ".", # score
                                            ".", # strand
                                            ".", # frame
                                            gtf_attributes # attribute
                                    ])

                                    continue
                            
                            # For novel transcripts
                            else:                         
                                if motif.location.start > cds_loc_position.end:
                                    break

                                # TODO check if conditional could be restricted to issubset() only, since it should include all the other conditions
                                # if motif.strand == transcript_exons_loc.strand and cds_positions[0] < motif.location.end and cds_positions[-1] > motif.location.start and motif_set.issubset(cds_positions):
                                if motif.strand == transcript_exons_loc.strand and transcript_exons_loc.start < motif.location.end and transcript_exons_loc.end > motif.location.start and motif_set.issubset(set(cds_loc_position)): # linea antigua
                                    # Until now all motif positions could be included inside the CDS, but the later could
                                    # contain a in-between fragment between them, so we perform a reverse search.
                                    filtered_cds = set(filter(lambda x: motif.location.start <= x < motif.location.end, sorted(cds_loc_position)))
                                    # Check that the motif is inside the transcript and has only one part (consecutive sequences)
                                    if filtered_cds.issubset(motif_set):
                                        try:    
                                            if motif.location.strand == 1:
                                                first_aa_position = int(cds_cm.g2p(motif.location.start))
                                            else:
                                            # For reverse strand,  take into account stop codon
                                                first_aa_position = int(cds_cm.g2p(motif.location.end-3))
                        
                                            last_aa_position = first_aa_position + len(motif_record.seq)
                                        except:
                                            continue

                                        # Extract the motif sequence and, if it matches the one in the input file, annotate it
                                        transcript_sequence = transcript_record.seq[first_aa_position:last_aa_position]
                                        if motif_record.seq == transcript_sequence:
                                            annotations = motif_record.annotations
                                            gtf_attributes = f'ID={annotations["id"]}; Name={annotations["name"]}; Desc={annotations["desc"]}_{annotations["database"]}'
                
                                            output_file_tsv.writerow([
                                                transcript_id, #seqname
                                                "UniProtKB/Swiss-Prot_Phosphosite", # source
                                                motif_record.annotations.get('feature'), # feature
                                                first_aa_position + 1, # start (1-based)
                                                last_aa_position, # end
                                                ".", # score
                                                ".", # strand
                                                ".", # frame
                                                gtf_attributes # attribute
                                        ])

                        # For ensembl/refseq transcriptomes, faster
                        else:
                            motif_protein_id = motif_record.name
                            # Keep version for EnsemblPlants
                            if "plants" in args.biomart_host and args.db.upper() == "ENSEMBL":
                                transcript_protein_id = transcript_record.name
                            else:
                                # Remove protein versions so they can match
                                transcript_protein_id = transcript_record.name.split(".")[0]

                            if motif_protein_id == transcript_protein_id:

                                try:    
                                    if motif.location.strand == 1:
                                        first_aa_position = int(cds_cm.g2p(motif.location.start))
                                    else:
                                        # For reverse strand,  take into account stop codon
                                        first_aa_position = int(cds_cm.g2p(motif.location.end-3))
                
                                    last_aa_position = first_aa_position + len(motif_record.seq)
                                except:
                                    continue

                                # Extract the motif sequence and, if it matches the one in the input file, annotate it
                                transcript_sequence = transcript_record.seq[first_aa_position:last_aa_position]
                                if motif_record.seq == transcript_sequence:
                                    annotations = motif_record.annotations
                                    gtf_attributes = f'ID={annotations["id"]}; Name={annotations["name"]}; Desc={annotations["desc"]}_{annotations["database"]}'

                                    output_file_tsv.writerow([
                                        transcript_id, #seqname
                                        "UniProtKB/Swiss-Prot_Phosphosite", # source
                                        motif_record.annotations.get('feature'), # feature
                                        first_aa_position + 1, # start (1-based)
                                        last_aa_position, # end
                                        ".", # score
                                        ".", # strand
                                        ".", # frame
                                        gtf_attributes # attribute
                                ])

    except Exception as ex:
        logging.error(str(ex))
        traceback.print_exc()
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
