def select_fasta_proteins(wildcards):
    return expand(rules.prepare_refseq_proteins.output, prefix=wildcards['prefix'])

def select_user_fasta_cdna(wildcards):
    return []
    
def select_fasta_cdna(wildcards):
    return expand(rules.prepare_refseq_cdna.output, prefix=wildcards['prefix'])[0]

def select_gtf(wildcards):
    return expand(rules.prepare_refseq_gtf.output.gtf, prefix=wildcards['prefix'])

def select_reference_gtf(wildcards):
    return expand(rules.prepare_refseq_gtf.output.gtf, prefix=wildcards['prefix'])

def select_sqanti_classification(wildcards):
    return expand(rules.run_refsqanti.output.classification, prefix=wildcards['prefix'], db=wildcards['db'])

def select_sqanti_output(wildcards):
    return rules.run_refsqanti.output

def select_nmd_file(wildcards):
    return expand(rules.run_refsqanti.output.NMD, prefix=wildcards['prefix'], db=wildcards['db'])

def select_prot_assoc(wildcards):
    return []

include: "Snakefile.smk"
