#!/usr/bin/env python3

'''
Author: Pedro Salguero - psalguero@cipf.es
Script to add a new field inside Description field that indicates the PosType for each annot (T, P, G or N)
'''

import sys

#line + (T,P,G,N)
def addPosType(res, line, posType):
    res.write(line[0:-1] + "; PosType=" + posType + "\n")


def updateGTF(filepath):
    # open new file
    res = open(filepath + "_mod", "w")
    # open annotation file and process all data
    with open(filepath, 'r') as f:
        # process all entries - no header line in file
        for line in f:
            if len(line) == 0:
                break
            else:
                if line and line[0] != "#":
                    fields = line.split("\t")
                    if len(fields) == 9:
                        if fields[1] == "tappAS":
                            if fields[2] == 'transcript':
                                addPosType(res, line, "T")
                            elif fields[2] == 'gene':
                                addPosType(res, line, "T")
                            elif fields[2] == 'CDS':
                                addPosType(res, line, "T")   
                            elif fields[2] == 'genomic':
                                addPosType(res, line, "G")
                            elif fields[2] == 'exon':
                                addPosType(res, line, "G")
                            elif fields[2] == 'splice_junction':
                                addPosType(res, line, "G")
                            elif fields[2] == 'protein':
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] == "COILS":   
                            if fields[2] == 'COILED':
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] == "GeneOntology":
                            if fields[2] in ('C', 'cellular_component'):
                                addPosType(res, line, "N")
                            elif fields[2] in ('F', 'molecular_function'):
                                addPosType(res, line, "N")
                            elif fields[2] in ('P', 'biological_process'):
                                addPosType(res, line, "N")
                            elif fields[2] in ('eco'):
                                addPosType(res, line, "N") # Fran tomato annot
                            else:
                                addPosType(res, line, "N")
                                #print(line)
                                #break

                        elif fields[1] == "MOBIDB_LITE": 
                            if fields[2] == 'DISORDER' :
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] == "NMD": 
                            if fields[2] == 'NMD':
                                addPosType(res, line, "T")
                            else:
                                print(line)
                                break

                        elif fields[1] in ("PAR-CLIP", "PAR-clip"):
                            if fields[2] in ('RNA_binding', 'RNA_Binding_Protein', 'RBP_Binding') or fields[2].startswith('RNA_binding_'):
                                addPosType(res, line, "T")
                            else:
                                print(line)
                                break

                        elif fields[1] == "PFAM":    
                            if fields[2] == 'DOMAIN':
                                addPosType(res, line, "P")
                            elif fields[2] in ("CLAN","clan"):
                                addPosType(res, line, "N")
                            else:
                                print(line)
                                break

                        elif fields[1] == "Provean": 
                            if fields[2] == 'FunctionalImpact':
                                addPosType(res, line, "N")
                            else:
                                print(line)
                                break

                        elif fields[1] in ("REACTOME","Reactome"):
                            if fields[2] in ('PATHWAY','pathway', 'Pathway'):
                                addPosType(res, line, "N")
                            else:
                                print(line)
                                break

                        elif fields[1] == "RepeatMasker":
                            if fields[2] == 'repeat':
                                addPosType(res, line, "T")
                            else:
                                print(line)
                                break

                        elif fields[1] == "SIGNALP_EUK": 
                            if fields[2] == 'SIGNAL':
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] == "TMHMM":   
                            if fields[2] == 'TRANSMEM':
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] == "TranscriptAttributes":
                            if fields[2] == '3UTR_Length':
                                addPosType(res, line, "T")
                            elif fields[2] == '5UTR_Length':
                                addPosType(res, line, "T")
                            elif fields[2] == 'CDS':
                                addPosType(res, line, "T")
                            elif fields[2] == 'polyA_Site':
                                addPosType(res, line, "T")
                            else:
                                print(line)
                                break

                        elif fields[1] == "UTRsite": 
                            if fields[2] == 'uORF':
                                addPosType(res, line, "T")
                            elif fields[2] == '5UTRmotif':
                                addPosType(res, line, "T")
                            elif fields[2] == 'PAS':
                                addPosType(res, line, "T")
                            elif fields[2] == '3UTRmotif':
                                addPosType(res, line, "T")
                            else:
                                print(line)
                                break

                        elif fields[1] in ("UniProtKB/Swiss-Prot_Phosphosite", "Swissprot_Phosphosite"):
                            if fields[2] == 'ACT_SITE':
                                addPosType(res, line, "P")
                            elif fields[2] == 'BINDING':
                                addPosType(res, line, "P")
                            elif fields[2] == 'PTM':
                                addPosType(res, line, "P")
                            elif fields[2] == 'MOTIF':
                                addPosType(res, line, "P")
                            elif fields[2] == 'COILED':
                                addPosType(res, line, "P")
                            elif fields[2] == 'TRANSMEM':
                                addPosType(res, line, "P")
                            elif fields[2] == 'COMPBIAS':
                                addPosType(res, line, "P")
                            elif fields[2] == 'INTRAMEM':
                                addPosType(res, line, "P")
                            elif fields[2] == 'NON_STD':
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] in ("cNLS_mapper", "NLS_mapper"): 
                            if fields[2] == 'MOTIF': 
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] in ("miRWalk", "mirWalk"): 
                            if fields[2] in('miRNA', 'miRNA_Binding'):
                                addPosType(res, line, "T")
                            else:
                                print(line)
                                break

                        elif fields[1] == "scanForMotifs":   
                            if fields[2] == 'PAS':
                                addPosType(res, line, "T")
                            elif fields[2] in ('3UTRmotif', "3'UTRmotif"):
                                addPosType(res, line, "T")
                            else:
                                print(line)
                                break

                        elif fields[1] == "MetaCyc":
                            if fields[2] == 'pathway':
                                addPosType(res, line, "N")
                            else:
                                print(line)
                                break

                        elif fields[1] == "KEGG":
                            if fields[2] in ('pathway','Pathway'):
                                addPosType(res, line, "N")
                            else:
                                print(line)
                                break

                        elif fields[1] == "SUPERFAMILY":
                            if fields[2] == 'DOMAIN':
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] == "SMART":
                            if fields[2] == 'DOMAIN':
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] == "TIGRFAM":
                            if fields[2] == 'DOMAIN':
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] == "psRNATarget":
                            if fields[2] == 'miRNA':
                                addPosType(res, line, "T")
                            else:
                                print(line)
                                break

                        elif fields[1] == "CORUM":
                            if fields[2] == 'Complex':
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] == "Orthologues":
                            if fields[2] == 'S.tuberosum':
                                addPosType(res, line, "N")
                            elif fields[2] in ('A.thaliana'):
                                addPosType(res, line, "N")
                            else:
                                print(line)
                                break

                        else:
                            print(line)
                            break

                    else:
                        raise Exception("Invalid line ({0}): {1}\n".format(idx, line))
                        break
        
        res.close()
        print("File has been created succesfuly")

def main():
    if len(sys.argv) > 2:
        sys.stderr.write("Usage: gen-matrix.py annotationFileName")
        raise SystemExit(1)
    if len(sys.argv) == 2:
        filepath = sys.argv[1]        
    print("Relabeling all features...")
    updateGTF(filepath)

if __name__ == "__main__":
    main()
