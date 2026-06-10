import sys

def select_fasta_proteins(wildcards):
    return rules.clean_sqanti_proteins.output

def select_user_cdna(wildcards):
    user_gtf = config.get("gtf_cdna", None)
    user_fasta = config.get("fasta_cdna", None)
    if not user_fasta and not user_gtf:
        sys.exit("FOR 'mytranscripts' DATABASE YOU NEED TO PROVIDE A CUSTOM FASTA WITH --config fasta_cdna=myfasta.fasta OR A CUSTOM GTF WITH --config gtf_cdna=mygtf.gtf")
    target = user_gtf if user_gtf else user_fasta
    try:
        import gzip
        opener = gzip.open if target.endswith(".gz") else open
        with opener(target, "rt") as f:
            first_line = ""
            for line in f:
                if line.strip():
                    first_line = line
                    break
            if user_fasta and not first_line.startswith(">"):
                sys.exit(f"ERROR: File '{target}' was provided as FASTA but doesn't start with '>'.\n FOR 'mytranscripts' DATABASE YOU NEED TO PROVIDE A CUSTOM FASTA WITH --config fasta_cdna=myfasta.fasta OR A CUSTOM GTF WITH --config gtf_cdna=mygtf.gtf")
            if user_gtf and first_line.startswith(">"):
                sys.exit(f"ERROR: File '{target}' was provided as GTF but looks like a FASTA (starts with '>').\n FOR 'mytranscripts' DATABASE YOU NEED TO PROVIDE A CUSTOM FASTA WITH --config fasta_cdna=myfasta.fasta OR A CUSTOM GTF WITH --config gtf_cdna=mygtf.gtf")
    except Exception as e:
        sys.exit(f"ERROR: Could not read file '{target}': {e}")
    return target

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
