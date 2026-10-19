#!/usr/bin/env python3
'''
Author: Lorena de la Fuente Lorente
Refactored: Carlos Martínez
Modified by Alessandra Martinez

Script that retrieves all the annotations produced by IsoAnnot
and produces a GFF3 file that can be read by tappAS.

'''

import argparse, sys, os, logging, csv, traceback, gc, re
import pandas as pd
from functools import reduce
from IsoAnnot import openfile, argparse_nullable, natural_key

# Order in which protein databases are used for the protein ID of a transcript (GitHub #10).
# SwissProt beats TrEMBL, which beats the native IDs. mytranscripts uses the Ensembl GTF as reference,
# so Ensembl comes before RefSeq there; the reference modes only list their own native database.
# "novel" is the Protein_<isoform> ID that transcript2reference.py gives to ORFs with no match.
PROTEIN_DB_RANK = {"swissprot": 0, "trembl": 1, "ensembl": 2, "refseq": 3, "novel": 4}


def resolve_protein_id(proteins):
    """
    Picks the protein ID of a transcript from the (protein_id, database) pairs that share its protein
    sequence, and returns it with the other IDs, which go to Desc. The result only depends on the set of
    pairs, so all transcripts that encode the same protein get the same ID and Desc in tappAS.
    Within a database, canonical UniProt entries come before isoform entries, even of a lower accession
    (P12345 before P12344-2), so a canonical entry is used whenever the sequence has one.
    """
    def _key(protein):
        protein_id, database = protein
        return (PROTEIN_DB_RANK.get(database, len(PROTEIN_DB_RANK)),
                re.search(r"-\d+$", protein_id) is not None, natural_key(protein_id))

    ranked = [protein_id for protein_id, _ in sorted(set(proteins), key=_key)]
    return ranked[0], ranked[1:]


def _free_df(df):
    del df
    gc.collect()


def get_structural_classification_df(classification_filename):
    '''
    Gets structural information for each transcript from the 
    SQANTI classification file.
    '''
    structural_classification = {}
    map_transcript_protein = {}

    with openfile(classification_filename, "rt") as read_handler:
        tsv_reader = csv.DictReader(read_handler, delimiter="\t")

        for line in tsv_reader:
            transcript_id = line.get("isoform").strip()
            transcript_location = None

            transcript_qualifiers = {
                "length": line.get("length"),
                "chrom": line.get("chrom"),
                "strand": line.get("strand"),
                "ref_transcript": line.get("associated_transcript").strip(),
                "ref_gene_name": line.get("associated_gene").replace(",", "_").strip(),
                "primary_class": line.get("structural_category").strip(),
                "second_class": ["RT_switching"] if line.get("RTS_stage") == "positive" else [],
                "cds_start": None if line.get("coding") != "coding" else line.get("CDS_start"),
                "cds_end": None if line.get("coding") != "coding" else line.get("CDS_end"),
                "ORF_length": None if line.get("coding") != "coding" else line.get("ORF_length"),
                "ref_protein": None,
                # "ref_CDS": None,
            }

            # Override CDS qualifiers for coding lines
            if line.get("coding") == "coding" and "refProt" in line:
                transcript_qualifiers["ref_protein"] = line.get("refProt")
                map_transcript_protein[line.get("refProt")] = transcript_id

            structural_classification[transcript_id] = transcript_qualifiers

    return structural_classification


def get_primary_dfs(classification_file, gene_desc_file, protein_assoc_info, extra_files):
    '''
    Gets gene, genomic (chr), transcript, CDS and protein information from
    SQANTI output into dataframes.
    '''

    # Get info from files
    transcript_collection = get_structural_classification_df(classification_file)
    # gene_descriptions = pd.read_csv(gene_desc_file, sep="\t", index_col=0, header=None) if gene_desc_file else None

    # Protein ID and Desc of each transcript, from the reference proteins with the same sequence
    # (GitHub #10). The file is empty when no protein matched.
    protein_ids = {}
    if protein_assoc_info and os.path.getsize(protein_assoc_info) > 0:
        protein_assoc = pd.read_csv(protein_assoc_info, sep="\t", index_col=None, header=None, dtype=str)
        protein_assoc.columns = ['transcript', 'protein', 'database', 'info']
        protein_assoc = protein_assoc[protein_assoc['info'] == 'sequence']
        for transcript_id, t_proteins in protein_assoc.groupby('transcript'):
            protein_ids[transcript_id] = resolve_protein_id(zip(t_proteins['protein'], t_proteins['database']))
    n_ref_protein_fallback = 0

    # Support info from GTF
    extra_second_files = None
    if extra_files is not None :
        extra_second_files = reduce(lambda left, right: pd.merge(left, right, on=0), [pd.read_csv(f, sep="\t", header=None, comment="#") for f in extra_files])
        extra_second_files.columns = ['transcript', 'info']

    transcript_df = [] 
    genomic_df = []
    protein_df = []

    # Retrieve information for each transcript in SQANTI classification file
    # Each row might require individual processing so we avoid vectorized operations.
    for transcript_id, t_prop in transcript_collection.items():
        transcript_chrom = t_prop.get("chrom")
        transcript_ref_id = transcript_id if t_prop.get("ref_transcript") == "novel" else t_prop.get("ref_transcript")
        transcript_attributes = [f"ID={transcript_ref_id}"]

        t_secondary_classes = t_prop.get("second_class")
        # t_gene_desc = gene_descriptions.at[transcript_id, 1]

        if extra_second_files is not None:
            t_extra_second = extra_second_files.loc[(extra_second_files['transcript'] == transcript_id) & (extra_second_files['info'].notna()), 'info'].tolist()
            t_secondary_classes.extend(t_extra_second)
            # TODO: iterate t_extra_second
            # for t_row in t_extra_second.iterrows():
            #     t_secondary_classes.extend(list(t_row)[1:])

        t_primary_class = t_prop.get("primary_class")
        if t_primary_class != "NA":
            transcript_attributes.append(f"primary_class={t_primary_class}")

        if len(t_secondary_classes):
            secondary_classes = ",".join(t_secondary_classes)
            transcript_attributes.append(f"secondary_class={secondary_classes}")

        # Append 'transcript' line
        transcript_df.append([transcript_id, "tappAS", "transcript", 1, t_prop.get("length"),
                          ".", t_prop.get("strand"), ".", "; ".join(transcript_attributes)])

        # Append 'gene' line
        t_gene_ref = t_prop.get("ref_gene_name")
        gene_desc = None #gene_descriptions.at[transcript_id, 1]

        gene_attributes = [f"ID={t_gene_ref}", f"Name={t_gene_ref}"]

        if gene_desc:
            gene_attributes.append(f"Desc={gene_desc}")

        transcript_df.append([transcript_id, "tappAS", "gene", 1, t_prop.get("length"),
                          ".", t_prop.get("strand"), ".", "; ".join(gene_attributes)])
        
        # Append 'genomic' line (Chromosome info)
        genomic_df.append([transcript_id, "tappAS", "genomic", 1, 1,
                           ".", t_prop.get("strand"), ".", f"Chr=chr{transcript_chrom}"])

        # Append 'CDS' line
        t_cds_start = t_prop.get("cds_start")

        if t_cds_start:
            if transcript_id in protein_ids:
                t_protein_id, t_other_proteins = protein_ids[transcript_id]
            elif t_prop.get("ref_protein") not in (None, "", "NA"):
                # Reference modes: refProt has no sequence in parsed_databases.json
                t_protein_id, t_other_proteins = t_prop["ref_protein"], []
                n_ref_protein_fallback += 1
            else:
                # IsoAnnotLite convention, instead of one "Novel" protein shared by all novel ORFs
                t_protein_id, t_other_proteins = "Protein_" + re.sub(r"[^A-Za-z0-9_.\-]", "_", transcript_id), []

            # Same attributes on the CDS and protein lines; no empty Desc (GitHub #10)
            cds_attributes = [f"ID={t_protein_id}"]
            if t_other_proteins:
                cds_attributes.append(f"Desc={','.join(t_other_proteins)}")

            transcript_df.append([transcript_id, "tappAS", "CDS", t_cds_start, t_prop.get("cds_end"),
                              ".", t_prop.get("strand"), ".", "; ".join(cds_attributes)])

            # Append 'protein' line. The end must be an integer: IsoAnnotLite crashes on '180.0' and
            # tappAS reads it as 0. Classifications made before this fix (GitHub #10) still have floats.
            protein_df.append([transcript_id, "tappAS", "protein", 1, int(float(t_prop.get("ORF_length"))),
                                  ".", t_prop.get("strand"), ".", "; ".join(cds_attributes)])

    if n_ref_protein_fallback:
        logging.warning(f"{n_ref_protein_fallback} coding transcripts have a refProt that is not in parsed_databases.json. "
                        f"Their refProt from the classification file is used as protein ID.")

    output_dfs = {
        "transcript": pd.DataFrame(transcript_df),
        "genomic": pd.DataFrame(genomic_df),
        "protein": pd.DataFrame(protein_df)
    }

    return output_dfs


def main():
    try:
        parser = argparse.ArgumentParser(description="IsoAnnot merging script")
        parser.add_argument('--output', required=True)
        parser.add_argument('--classification_file', required=True)
        parser.add_argument('--extra_files', nargs='*', required=False, type=argparse_nullable)
        parser.add_argument('--protein_association', nargs='?', default=False, required=False, type=argparse_nullable)
        parser.add_argument('--gene_desc_file', required=False, nargs='?', const=None, type=argparse_nullable)
        parser.add_argument('--input_transcripts', nargs='+', required=True)
        parser.add_argument('--input_genomic', nargs='+', required=False)
        parser.add_argument('--input_protein', nargs='+', required=False)

        args = parser.parse_args()

        # Primary data frames
        logging.info("Generating basic GTF structure")
        block_dfs = get_primary_dfs(args.classification_file, args.gene_desc_file, args.protein_association, args.extra_files)

        def _read_gtf_files(file_list, filter_ids):
            gtf_dfs = []
            for gtf_file in file_list:
                try:
                    logging.info(f"Reading GTF file {gtf_file}...")
                    gtf_contents = pd.read_csv(gtf_file, sep="\t", header=None)
                    n_lines = len(gtf_contents)

                    gtf_contents = gtf_contents[gtf_contents[0].isin(filter_ids)]
                    # A version mismatch used to empty a layer without any message
                    if len(gtf_contents) < n_lines:
                        logging.warning(f"{n_lines - len(gtf_contents)} of {n_lines} lines in {gtf_file} were dropped because "
                                        f"their transcript ID is not in the classification file. If that is most of "
                                        f"the file, check the transcript_versioned setting of the config file.")

                    gtf_dfs.append(gtf_contents)
                except Exception as e:
                    # Ignore empty file errors
                    logging.error(f"Error reading GTF file {gtf_file}: " + str(e) + ". Ignoring.")
                    continue

            # pd.concat([]) fails when every file of a block is empty
            return gtf_dfs if gtf_dfs else [pd.DataFrame()]

        # Get transcript IDs annotated in primary DFs
        existing_ids = block_dfs.get("transcript")[block_dfs.get("transcript").columns[0]]

        # Read Transcript GTFs
        logging.info("Reading transcript block GTF files")
        transcript_dfs = pd.concat(_read_gtf_files(args.input_transcripts, existing_ids), sort=False, ignore_index=True)

        # Read Genomic GTFs
        logging.info("Reading genomic block GTF files")
        genomic_dfs = pd.concat(_read_gtf_files(args.input_genomic, existing_ids), sort=False, ignore_index=True)

        # Read Protein GTFs
        logging.info("Reading protein block GTF files")
        protein_dfs = pd.concat(_read_gtf_files(args.input_protein, existing_ids), sort=False, ignore_index=True)

        # Create the final GTF using the required tappAS structure: 3 blocks
        logging.info("Starting merging process")

        final_annotation = block_dfs.get("transcript")

        # Append the "transcript" information
        logging.info("Adding transcript block...")
        final_annotation = pd.concat([final_annotation, transcript_dfs], sort=False, ignore_index=True)

        # Free old dfs
        _free_df(transcript_dfs)

        # Append the "genomic" block
        logging.info("Adding genomic block...")
        final_annotation = pd.concat([final_annotation, block_dfs.get("genomic"), genomic_dfs], sort=False, ignore_index=True)

        # Free old dfs
        _free_df(genomic_dfs)

        # Append the "protein" block
        logging.info("Adding protein block...")
        final_annotation = pd.concat([final_annotation, block_dfs.get("protein"), protein_dfs], sort=False, ignore_index=True)

        # Free old dfs
        _free_df(protein_dfs)
        _free_df(block_dfs)

        # Sort by first column (transcript id)
        logging.info("Sorting by transcript...")

        final_annotation.rename(columns={final_annotation.columns[0]: "transcript"}, inplace = True)
        final_annotation = final_annotation.rename_axis('myindex').sort_values(by= ["transcript", "myindex"], ascending = [True, True], axis=0)

        logging.info("Writing final GTF file...")
        final_annotation.to_csv(args.output, sep="\t", header=False, index=False)

    except Exception as ex:
        logging.error(str(ex))
        traceback.print_exc()
        sys.exit(os.EX_SOFTWARE)

    sys.exit(os.EX_OK)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
