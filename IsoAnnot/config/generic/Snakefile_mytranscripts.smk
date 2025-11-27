def select_fasta_proteins(wildcards):
    return rules.clean_sqanti_proteins.output

def select_user_fasta_cdna(wildcards):
    user_fasta = config.get("fasta_cdna", None)
    if not user_fasta:
        print("FOR 'mytranscripts' DATABASE YOU NEED TO PROVIDE A CUSTOM FASTA WITH --config fasta_cdna=myfasta.fasta")
        sys.exit(os.EX_SOFTWARE)
    return user_fasta

def select_fasta_cdna(wildcards):
    return rules.run_sqanti.output.corrected_cdna

def select_gtf(wildcards):
    return rules.run_sqanti.output.gtf

def select_reference_gtf(wildcards):
    return rules.prepare_ensembl_gtf.output

def select_sqanti_classification(wildcards):
    return rules.run_sqanti.output.classification

def select_sqanti_output(wildcards):
    return rules.run_sqanti.output

def select_nmd_file(wildcards):
    return rules.transcript_to_reference.output.nmd

def select_prot_assoc(wildcards):
    return rules.transcript_to_reference.output.protein_assoc

include: "Snakefile.smk"
