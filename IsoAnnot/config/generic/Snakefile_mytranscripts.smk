import sys, os

def select_fasta_proteins(wildcards):
    return "data/Hsapiens/output/mytranscripts/sqanti_fasta_proteins.fasta"

def select_user_fasta_cdna(wildcards):
    user_fasta = config.get("fasta_cdna", None)
    if not user_fasta:
        print("FOR 'mytranscripts' DATABASE YOU NEED TO PROVIDE A CUSTOM FASTA WITH --config fasta_cdna=myfasta.fasta")
        sys.exit(os.EX_SOFTWARE)
    return user_fasta

def select_fasta_cdna(wildcards):
    return "data/Hsapiens/output/mytranscripts/sqanti_corrected.fasta"

def select_gtf(wildcards):
    return "data/Hsapiens/output/mytranscripts/sqanti_gtf.gtf"

def select_reference_gtf(wildcards):
    return rules.prepare_ensembl_gtf.output

def select_sqanti_classification(wildcards):
    return "data/Hsapiens/output/mytranscripts/sqanti_classification.txt"

def select_sqanti_output(wildcards):
    return rules.run_sqanti.output

def select_nmd_file(wildcards):
    return expand(rules.transcript_to_reference.output.nmd, db=wildcards['db'])

def select_prot_assoc(wildcards):
    return rules.transcript_to_reference.output.protein_assoc

include: "Snakefile.smk"
