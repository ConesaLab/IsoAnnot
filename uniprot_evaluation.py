#!/usr/bin/env python3

'''
Author: Alessandra Martínez
Script to compare the protein annotations produced by IsoAnnot to those from UniprotKB.
'''

import gff3_parser, pandas

input_uniprot = "drerio/UP000000437_7955.dat"
input_uniprot_parsed = "drerio/uniprot_parsed.txt"
input_isoannot_gff3 = "drerio/layer_uniprot.gtf"
valid_features = ["ACT_SITE", "BINDING", "CA_BIND", "DNA_BIND", "MOD_RES", "CARBOHYD", "LIPID",
                        "DISULFID", "METAL", "NONT_STD", "COILED", "INTRAMEM", "MOTIF", "REGION", "SITE",
                        "COMPBIAS", "CROSSLNK", "TRANSMEM", "ZN_FING"]

# Get uniprot canonical isoforms (look for Sequence=Displayed; in uniprot .dat file)
uniprot_principal_isoforms = []
uniprot_principal_accessions = []

with open(input_uniprot, "r") as uniprot_input:
    handler1 = uniprot_input.readlines()
    for line in handler1:
        if "Sequence=Displayed;" in line:
            uniprot_principal_isoforms.append(line.split(";")[0].split("=")[1])
            uniprot_principal_accessions.append(line.split(";")[0].split("=")[1].split("-")[0])

print("Uniprot principal protein isoforms: " + str(len(uniprot_principal_accessions)))

print(uniprot_principal_isoforms)
# Get transcript identifiers (ENST) associated to those uniprot isoforms
transcript_principal_isoforms = []
features_to_eval = []
transcript_protein_dict = {}

with open(input_uniprot_parsed, "r") as uniprot_parsed_input:
    handler2 = uniprot_parsed_input.readlines()
    for line in handler2:
        # if line.split("\t")[2] == "DR" and line.split("\t")[3] == "EnsemblPlants":
        # if line.split("\t")[2] == "DR" and line.split("\t")[3] == "EnsemblMetazoa": 
        if line.split("\t")[2] == "DR" and line.split("\t")[3] == "Ensembl": 
            if "[" in line and line.split("[")[1].split("]")[0] in uniprot_principal_isoforms:
                # transcript_principal_isoforms.append(line.split("\t")[4])
                transcript_principal_isoforms.append(line.split("\t")[4]).split(".")[0]  # Remove version
                # transcript_protein_dict[line.split("\t")[4]] = line.split("[")[1].split("]")[0]
                transcript_protein_dict[line.split("\t")[4].split(".")[0]] = line.split("[")[1].split("]")[0] # Remove version

        # Get features annotated by isoannot for those uniprot accessions
        if line.split("\t")[2] == "FT":
            if line.split("\t")[3] in valid_features and line.split("\t")[0] in uniprot_principal_accessions:
                features_to_eval.append(line)

print("Ensembl principal isoforms: " + str(len(transcript_principal_isoforms)))

features_to_eval_sep = []
for line in features_to_eval:
    line = line.rstrip().replace(".", "\t").split("\t")
    features_to_eval_sep.append(line[0:7])

features_to_eval_df = pandas.DataFrame(features_to_eval_sep, columns= ["Uniprot", "Species", "FT", "Type", "Start", "End", "Desc"])

#Filter uniprot annotations and principal isoforms transcripts from isoAnnot GFF3
isoannot_df = gff3_parser.parse_gff3(input_isoannot_gff3,parse_attributes=True)
isoannot_df = isoannot_df[isoannot_df['Source']=="UniProtKB/Swiss-Prot_Phosphosite"]
isoannot_df = isoannot_df[isoannot_df[' Desc'].str.contains("SwissProt")]
isoannot_df = isoannot_df[isoannot_df['Seqid'].isin(transcript_principal_isoforms)]

no_annot = 0
proteins_not_annotated = []
fewer_annot = 0
proteins_fewer = []
more_annot = 0
proteins_more = []
wrong_positions = 0
proteins_wrong = []
wrong_protein = 0
ok_protein = 0

# For each pair uniprot_protein:ensembl_transcript compare the number of features annotated
# by uniprot and isoannot, and the positions
for element in transcript_protein_dict.items():
    transcript = element[0]
    protein = element[1].split("-")[0]
    original_annotation = features_to_eval_df[features_to_eval_df["Uniprot"]==protein]
    # Uniprot has duplicated entries (eg: Q6NNF2 binding site in 683)
    original_annotation = original_annotation[['Type','Start', 'End', 'Desc']].drop_duplicates()
    isoannot_annotation = isoannot_df[isoannot_df["Seqid"]==transcript] 
    isoannot_annotation = isoannot_annotation[['Type','Start', 'End',  ' Desc']].drop_duplicates()

    if original_annotation.shape[0] < isoannot_annotation.shape[0]:
        more_annot +=1
        proteins_more.append(transcript)
    elif original_annotation.shape[0] > isoannot_annotation.shape[0]:
        if isoannot_annotation.shape[0] == 0:
            no_annot +=1
            proteins_not_annotated.append(transcript)
        else:
            fewer_annot +=1
            proteins_fewer.append(transcript)

    # If same number of features, compare positions
    else:
        # Uniprot is 0 based, isoannot is 1 based
        original_start = original_annotation["Start"].astype(int) +1
        original_start = original_start.sort_values().to_list()
        original_end = original_annotation["End"].astype(int).sort_values().to_list()
        isoannot_start = isoannot_annotation["Start"].astype(int).sort_values().to_list()
        isoannot_end = isoannot_annotation["End"].astype(int).sort_values().to_list()

        if isoannot_start == original_start and isoannot_end == original_end:
            ok_protein += 1
        else:
            wrong_positions +=1
            proteins_wrong.append(transcript)

print("Isoforms not annotated: " + str(no_annot))
print("Isoforms with less annotations than Uniprot: " + str(fewer_annot))
print("Isoforms with more annotations than Uniprot: " + str(more_annot))
print("Isoforms with same number of annotations but incorrect positions: " + str(wrong_positions))
print("Correctly annotated isoforms: " + str(ok_protein))

# print("Proteins with no annotations")
# print(proteins_not_annotated[0:10])
# print("Proteins with fewer annotations")
# print(proteins_fewer[0])
# print("Proteins with more annotations")
# print(proteins_more[0])
# print("Proteins with wrong positions")
# print(proteins_wrong[0])
