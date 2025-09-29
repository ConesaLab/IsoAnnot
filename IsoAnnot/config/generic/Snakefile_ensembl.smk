def select_fasta_proteins(wildcards):
    return rules.get_ensembl_proteins.output.fa

def select_user_fasta_cdna(wildcards):
    return []

def select_fasta_cdna(wildcards):
    return rules.prepare_ensembl_cdna.output.fa

def select_gtf(wildcards):
    return rules.prepare_ensembl_gtf.output

def select_reference_gtf(wildcards):
    return rules.prepare_ensembl_gtf.output

# def select_extra_files(wildcards):
#     return expand(rules.prepare_ensembl_support.output.transcripts, prefix=wildcards['prefix'])

def select_sqanti_classification(wildcards):
    return expand(rules.run_refsqanti.output.classification, db=wildcards['db'])

def select_sqanti_output(wildcards):
    return rules.run_refsqanti.output

def select_nmd_file(wildcards):
    return expand(rules.run_refsqanti.output.NMD, db=wildcards['db'])

def select_prot_assoc(wildcards):
    return []

include: "Snakefile.smk"
