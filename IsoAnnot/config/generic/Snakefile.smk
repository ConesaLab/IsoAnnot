import os

configfile: "config/generic/config.yaml"
prefix = config["prefix"]
db = config["db"]
species_name = config["species_name"]

def _output_layer_db(layer_name, external_rule=[], wildcards=None):

    if callable(external_rule):
        external_rule = external_rule({"prefix": prefix, "db": db })
    if len(config.get(layer_name, [])) or len(external_rule):
        return os.path.join("data", prefix, "output", "{db}", "layers", f"{layer_name}.gtf")
    else:
        return []

def _optional_param(param_name, param_option):
    config_value = config.get(param_name, None)

    if config_value:
        output_param = f"{param_option} {config_value}"
    else:
        output_param = ""

    return output_param

def _remove_extension(config_value, all=False):
    filename = os.path.basename(config_value)

    if config_value.endswith("gz"):
        output_name = filename[:-7] if all else filename[:-3]
    elif all:
        ouput_name = filename[:-4]

    return output_name


rule all:
    input:
        expand("data/{prefix}/{species_name}_tappas_{db}_annotation_file.gff3_mod",
               prefix=prefix,
               species_name=species_name,
               db=db)

# GET AND PREPARE

# rule get_conversion_file:
#     output:
#         expand("data/{{prefix}}/config/ensembl/{filename}", filename = os.path.basename(config["conversion_file"]))[0]
#     params:
#         URL=config["conversion_file"]

#     shell:
#         """
#         wget -nv -P data/{wildcards.prefix}/config/ensembl/ {params.URL}
#         """

# rule prepare_conversion_file:
#     conda:
#         "../../envs/isoannotpy.yaml"
#     input:
#         rules.get_conversion_file.output
#     output:
#         expand("data/{{prefix}}/config/ensembl/{filename}", filename = _remove_extension(config["conversion_file"]))[0]
#     shell:
#         """
#         gunzip -k {input}
#         """

rule get_refseq_acc: 
    output: 
        os.path.join("data",prefix,"config","refseq", os.path.basename(config["refseq_chr_accessions"])) 
    params: 
        URL=config["refseq_chr_accessions"] 
    log: 
        os.path.join("logs", prefix, "get_refseq_acc.log")
    shell: 
        """ 
        wget -nv -P data/{prefix}/config/refseq/ {params.URL} &> {log} 
        """


checkpoint get_refseq_proteins:
    output:
        directory(os.path.join("data", prefix, "config", "refseq", "faa"))
    params:
        URL=config["refseq_protein_dir"]
    log:
        os.path.join("logs", prefix, "get_refseq_proteins.log")
    shell:
        """
        wget -nv -P data/{prefix}/config/refseq/faa/ -nd -np -r -nH -l1 -A *protein.faa.gz {params.URL} &> {log}
        """

def gather_refseq_proteins(wildcards):
    checkpoint_output = checkpoints.get_refseq_proteins.get(**wildcards).output[0]
    fname_vals = glob_wildcards(os.path.join(checkpoint_output, "{filename}.faa.gz")).filename
    return expand(os.path.join(checkpoint_output, "{filename}.faa.gz"), dir=checkpoint_output, filename=fname_vals)[0]


rule prepare_refseq_proteins:
    input:
       gather_refseq_proteins
    output:
        os.path.join("data", prefix, "config", "refseq", "all_refseq_proteins.faa")
    log:
        os.path.join("logs", prefix, "prepare_refseq_proteins.log")
    shell:
       """
       zcat data/{prefix}/config/refseq/faa/*.faa.gz > {output} 2> {log}
       """

checkpoint get_refseq_cdna:
    output:
        directory(os.path.join("data", prefix, "config", "refseq", "fna"))
    params:
        URL=config["refseq_protein_dir"]
    log:
        os.path.join("logs", prefix, "get_refseq_cdna.log")
    shell:
        """
        wget -nv -P data/{prefix}/config/refseq/fna/ -nd -np -r -nH -l1 -A *rna.fna.gz {params.URL} &> {log}
        """


def gather_refseq_cdna(wildcards):
    checkpoint_output = checkpoints.get_refseq_cdna.get(**wildcards).output[0]
    fname_vals = glob_wildcards(os.path.join(checkpoint_output, "{filename}.fna.gz")).filename
    return expand(os.path.join(checkpoint_output, "{filename}.fna.gz"), dir=checkpoint_output, filename=fname_vals)[0]


rule prepare_refseq_cdna:
    input:
        gather_refseq_cdna
    output:
        fa=os.path.join("data", prefix, "config", "refseq", "all_refseq_cdna.fna")
    log:
        os.path.join("logs", prefix, "prepare_refseq_cdna.log")
    shell:
        """
        zcat data/{prefix}/config/refseq/fna/*.fna.gz > {output} 2> {log}
        """

rule get_refseq_gtf:
    output:
        os.path.join("data", prefix, "config", "refseq", os.path.basename(config["refseq_gtf"]))
    params:
        URL=config["refseq_gtf"]
    log:
        os.path.join("logs", prefix, "get_refseq_gtf.log")
    shell:
        """
        wget -nv -P data/{prefix}/config/refseq/ {params.URL} &> {log}
        """


rule prepare_refseq_gtf:
    input:
        rules.get_refseq_gtf.output
    output:
        os.path.join("data", prefix, "config", "refseq", f"{_remove_extension(config['refseq_gtf'], all=True)}_nopartial_nc.gtf")
    log:
        os.path.join("logs", prefix, "prepare_refseq_gtf.log")
    shell:
        """
        zcat {input}|grep -v -F 'partial=true' | grep NC_0  > {output} 2> {log}
        """


rule get_ensembl_proteins:
    conda:
        "../../envs/isoannotpy.yaml"   
    output:
        gz=os.path.join("data", prefix, "config", "ensembl", os.path.basename(config["ensembl_proteins"])),
        fa=os.path.join("data", prefix, "config", "ensembl", _remove_extension(config["ensembl_proteins"]))
    params:
        URL=config["ensembl_proteins"]
    log:
        os.path.join("logs", prefix, "get_ensembl_proteins.log")
    shell:
        """
        (wget -nv -P data/{prefix}/config/ensembl/ {params.URL}
        gunzip -c {output.gz} > {output.fa}.tmp
        sed s/\*//g {output.fa}.tmp > {output.fa}
        rm {output.fa}.tmp) 2> {log}
        """
     

rule get_ensembl_cdna:
    output:
        os.path.join("data", prefix, "config", "ensembl", os.path.basename(config["ensembl_cdna"]))
    params:
        URL=config["ensembl_cdna"]
    log:
        os.path.join("logs", prefix, "get_ensembl_cdna.log")
    shell:
        """
        wget -nv -P data/{prefix}/config/ensembl/ {params.URL} &> {log}
	"""


rule prepare_ensembl_cdna:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.get_ensembl_cdna.output
    output:
        fa=os.path.join("data", prefix, "config", "ensembl", _remove_extension(config["ensembl_cdna"]))
    log:
        os.path.join("logs", prefix, "prepare_ensembl_cdna.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """


rule get_ensembl_reference:
    output:
        os.path.join("data", prefix, "config", "ensembl", os.path.basename(config["ensembl_reference"]))
    params:
        URL=config["ensembl_reference"]
    log:
        os.path.join("logs", prefix, "get_ensembl_reference.log")
    shell:
        """
        wget -nv -P data/{prefix}/config/ensembl/ {params.URL} &> {log}
        """

rule prepare_ensembl_reference:
    input:
        rules.get_ensembl_reference.output
    output:
        os.path.join("data", prefix, "config", "ensembl", _remove_extension(config["ensembl_reference"]))
    log:
        os.path.join("logs", prefix, "prepare_ensembl_reference.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """

rule get_ensembl_gtf:
    output:
        os.path.join("data", prefix, "config", "ensembl", os.path.basename(config["ensembl_gtf"]))
    params:
        URL=config["ensembl_gtf"]
    log:
        os.path.join("logs", prefix, "get_ensembl_gtf.log")
    shell:
        """
        wget -nv -P data/{prefix}/config/ensembl/ {params.URL} &> {log}
        """


rule prepare_ensembl_gtf:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.get_ensembl_gtf.output
    output:
        os.path.join("data", prefix, "config", "ensembl", _remove_extension(config["ensembl_gtf"]))
    log:
        os.path.join("logs", prefix, "prepare_ensembl_gtf.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """


rule get_pfam_clan:
    output:
        os.path.join("data", "global", "pfam", os.path.basename(config["pfam_clan_url"]))
    params:
        URL=config["pfam_clan_url"]
    log:
        os.path.join("logs", prefix, "get_pfam_clan.log")
    shell:
        """
        wget -nv -P data/global/pfam/ {params.URL} &> {log}
        """


rule prepare_pfam_clan: 
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.get_pfam_clan.output
    output:
        os.path.join("data", "global", "pfam", _remove_extension(config["pfam_clan_url"]))
    log:
        os.path.join("logs", prefix, "prepare_pfam_clan.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """


rule get_uniprot_data: 
    output:
        [os.path.join("data",prefix,"config","uniprot",os.path.basename(uniprot_url)) for uniprot_url in config["uniprot_dat"] + config["uniprot_fasta"]]
    log:
        [os.path.join("logs", prefix, "get_uniprot_data", f"{_remove_extension(uniprot_url)}.log") for uniprot_url in config["uniprot_dat"] + config["uniprot_fasta"]]
    run:
        for uniprot_url, log_file in zip(config["uniprot_dat"] + config["uniprot_fasta"], log):
            shell(f"wget -nv -P data/{prefix}/config/uniprot/ {uniprot_url} &> {log_file}")


rule prepare_uniprot_data:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        [file for file in rules.get_uniprot_data.output if file.endswith(".dat.gz")]
    output:
        os.path.join("data", prefix, "config", "uniprot", "uniprot_parsed.txt")
    log:
        os.path.join("logs", prefix, "prepare_uniprot_data.log")
    shell:
        """
        scripts/uniprot_parse.py --uniprot_files {input} --output {output} &> {log}
        """


rule get_reactome:
    output:
        os.path.join("data", "global", os.path.basename(config["reactome"]))
    params:
        URL=config["reactome"]
    log:
        os.path.join("logs", prefix, "get_reactome.log")
    shell:
        """
        wget -nv --no-check-certificate -P data/global/ {params.URL} &> {log}
        """
	

# RUN

rule run_refsqanti:
    conda:
        "../../envs/sqanti.yaml"
    input:
        gtf=select_reference_gtf,
        reference=rules.prepare_ensembl_reference.output,
        chr_ref=rules.get_refseq_acc.output  # we use this even in ensembl mode to filter only chromosome sequences (avoid MT)
    output:
        classification=os.path.join("data", prefix, "config", "{db}", "sqanti_classification.txt"),
        junctions=os.path.join("data", prefix, "config", "{db}", "sqanti_junctions.txt"),
        NMD=os.path.join("data", prefix, "config", "{db}", "sqanti_NMD.txt")
    log:
        os.path.join("logs", prefix, "{db}", "run_refsqanti.log")
    shell:
        """
        scripts/referenceSQANTI.py --gtf_file {input.gtf} --reference_file {input.reference} --chr_ref {input.chr_ref} --database {wildcards.db} --output_classification {output.classification} --output_junctions {output.junctions} --output_nmd {output.NMD} &> {log}
        """

rule run_gmap_index:
    conda:
        "../../envs/sqanti.yaml"
    input:
        ref_genome=rules.prepare_ensembl_reference.output
    output:
        directory(os.path.join("data", prefix, "config", "ensembl", "gmap_index"))
    params:
        outdir=os.path.join("data",prefix,"config","ensembl"),
        index_name="gmap_index"
    log:
        os.path.join("logs", prefix, "run_gmap_index.log")
    shell:
        """
        gmap_build -D {params.outdir} -d {params.index_name} {input.ref_genome} &> {log}
        """


rule run_sqanti:
    conda:
        "../../envs/sqanti.yaml"
    input:
        user_fasta=select_user_fasta_cdna,
        reference_gtf=select_reference_gtf,
        genome_fasta=rules.prepare_ensembl_reference.output,
        genome_fasta_index=rules.run_gmap_index.output
    output:
        corrected_cdna=os.path.join("data", prefix, "output", "{db}", "sqanti_corrected.fasta"),
        fasta_proteins=os.path.join("data", prefix, "output", "{db}", "sqanti_fasta_proteins.fasta"),
        gtf=os.path.join("data", prefix, "output", "{db}", "sqanti_gtf.gtf"),
        classification=os.path.join("data", prefix, "output", "{db}", "sqanti_classification.txt"),
        junctions=os.path.join("data", prefix, "output", "{db}", "sqanti_junctions.txt"),
    params:
        outdir=os.path.join("data",prefix,"output","{db}"),
        out_name="sqanti"
    log:
        os.path.join("logs", prefix, "{db}", "run_sqanti.log")
    shell:
        """
        scripts/sqanti/sqanti_qc.py {input.user_fasta} {input.reference_gtf} {input.genome_fasta} -d {params.outdir} -x {input.genome_fasta_index}/gmap_index -o {params.out_name} -n 8 --skipORF &> {log}
        """


rule run_utrscan:
    input:
        select_fasta_cdna
    output:
        os.path.join("data", prefix, "output", "{db}", "utrscan.txt")
    log:
        os.path.join("logs", prefix, "{db}", "run_utrscan.log")
    shell:
        """
        software/bin/UtrScan -SIGNALLIST -COMMAND=software/bin/UtrSite.Command -INPUT={input} -OUTPUT={output} &> {log}
        """

rule run_repeatmasker: 
    conda:
        "../../envs/repeats.yaml"
    input:
        select_fasta_cdna
    output:
        os.path.join("data", prefix, "output", "{db}", "repeat_masker", f"{os.path.basename(select_fasta_cdna(config))}.out")
    params:
        species_name=config["species_name"],
        outdir=os.path.join("data", prefix, "output", "{db}", "repeat_masker/")
    log:
        os.path.join("logs", prefix, "{db}", "run_repeatmasker.log")
    shell:
        r"""
        LIBDIR="$CONDA_PREFIX/lib"
        if [ ! -e "$LIBDIR/libnsl.so.1" ]; then
            ln -s "$LIBDIR/libnsl.so.3" "$LIBDIR/libnsl.so.1"  #se crea un enlace simbólico, no hay ninguna versión que instale la dependencia requerida
        fi
        export LD_LIBRARY_PATH="$LIBDIR:$LD_LIBRARY_PATH"
        RepeatMasker {input} -species {params.species_name} -dir {params.outdir} &> {log}
        """

 
checkpoint run_interproscan:
    conda:
        "../../envs/interpro_java.yaml"
    input:
        select_fasta_proteins
    output:
        directory(os.path.join("data", prefix, "output", "{db}", "interproscan", "splitProteins"))  
    params:
        interproscan_path=config["interproscan_path"]
    log:
        os.path.join("logs", prefix, "{db}", "run_interproscan.log")
    shell:
        """
        mkdir -p {output} && {params.interproscan_path} -i {input} -d {output} --disable-precalc  -appl Coils,Pfam,MobiDBLite,SignalP_EUK,TMHMM  -f XML -iprlookup &> {log}
        """

def gather_interproscan(wildcards):
    checkpoint_output=checkpoints.run_interproscan.get(**wildcards).output[0]
    fname_vals = glob_wildcards(os.path.join(checkpoint_output, "{filename}.xml")).filename
    return expand(os.path.join(checkpoint_output, "{filename}.xml"), dir=checkpoint_output, filename=fname_vals)[0]


rule parse_interproscan:
    conda:
        "../../envs/sqanti.yaml"
    input:
        gather_interproscan
    output:
        os.path.join("data", prefix, "output", "{db}", "interproscan", "interproResultsPfam.tsv")
    log:
        os.path.join("logs", prefix, "{db}", "parse_interproscan.log")
    shell:
        """
        scripts/parseInterproscanXml.py --interproscan_files {input} --output {output} &> {log}
        """


rule run_gtftogenepred:
    conda:
        "../../envs/genepred.yaml"
    input:
        select_gtf
    output:
        os.path.join("data", prefix, "config", "{db}", "genePrediction.txt")
    log:
        os.path.join("logs", prefix, "{db}", "run_gtftogenepred.log")
    shell:
        """
        gtfToGenePred {input} {output} -genePredExt -allErrors -ignoreGroupsWithoutExons &> {log}
        """


rule get_genomic_coordinates:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        uniprot_fasta = [file for file in rules.get_uniprot_data.output if file.endswith(".fasta.gz")],
        refseq_fasta = rules.prepare_refseq_proteins.output if config["refseq_protein_fasta"] else [],
        ensembl_fasta = rules.get_ensembl_proteins.output.fa,
        phosphosite_files = [
            "data/global/PSP_data/Acetylation_site_dataset",
            "data/global/PSP_data/Methylation_site_dataset",
            "data/global/PSP_data/O-GalNAc_site_dataset",
            "data/global/PSP_data/O-GlcNAc_site_dataset",
            "data/global/PSP_data/Phosphorylation_site_dataset",
            "data/global/PSP_data/Sumoylation_site_dataset",
            "data/global/PSP_data/Ubiquitination_site_dataset"],
        refseq_gtf = rules.prepare_refseq_gtf.output if config["refseq_gtf"] else [],
        ensembl_gtf = rules.prepare_ensembl_gtf.output,
        uniprot_parsed = rules.prepare_uniprot_data.output,
        chr_ref = rules.get_refseq_acc.output
    output:
        protein=os.path.join("data", prefix, "config", "uniprot", "uniprot_gcord_proteinGenomic.txt"),
        domain=os.path.join("data", prefix, "config", "uniprot", "uniprot_gcord_domainGenomic.txt"),
    params:
        biomart_host = config.get("biomart_host", [])
    log:
        os.path.join("logs", prefix, "get_genomic_coordinates.log")
    shell:
        """
        scripts/uniprotPhosphosite_genomicCoordinates.py --uniprot_fasta {input.uniprot_fasta} --refseq_fasta {input.refseq_fasta} --ensembl_fasta {input.ensembl_fasta} --uniprot_parsed {input.uniprot_parsed} --phosphosite_files {input.phosphosite_files} --refseq_gtf {input.refseq_gtf} --ensembl_gtf {input.ensembl_gtf}  --output_protein {output.protein} --output_domain {output.domain} --chr_ref {input.chr_ref} --biomart_host {params.biomart_host} &> {log}
        """
        

rule get_uniprot_phosphosite_annotation:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        fasta_orf=select_fasta_proteins,
        classification=select_sqanti_classification,
        gene_prediction=rules.run_gtftogenepred.output,
        motif_info=rules.get_genomic_coordinates.output.domain,
        chr_ref = rules.get_refseq_acc.output,
        protein_assoc=select_prot_assoc
    params:
        keep_version = config["transcript_versioned"],
        biomart_host = config.get("biomart_host", [])        
    output:
        os.path.join("data", prefix, "output", "{db}", "uniprot_Phosphosite_info.txt")
    log:
        os.path.join("logs", prefix, "{db}", "get_uniprot_phosphosite_annotation.log")
    shell:           
        """
        scripts/uniprotPhosphosite_annotation.py --orf_fasta {input.fasta_orf} --classification_file {input.classification} --genepred_file {input.gene_prediction} --uniprotmotif_file {input.motif_info} --chr_ref {input.chr_ref} --protein_association {input.protein_assoc} --keep_version {params.keep_version} --db {wildcards.db} --biomart_host {params.biomart_host} --output {output} &> {log}
        """


rule parse_protein_databases:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        uniprot_fasta = [file for file in rules.get_uniprot_data.output if file.endswith(".fasta.gz")],
        refseq_fasta=rules.prepare_refseq_proteins.output,
        ensembl_fasta=rules.get_ensembl_proteins.output.fa,
    output:
        os.path.join("data", prefix, "config", "parsed_databases.json")
    params:
        ensembl_fasta_regex=_optional_param("ensembl_fasta_regex", "--ensembl_fasta_regex"),
        refseq_fasta_regex=_optional_param("refseq_fasta_regex", "--refseq_fasta_regex"),
    log:
        os.path.join("logs", prefix, "parse_protein_databases.log")
    shell:
        """
        scripts/parse_protein_databases.py --uniprot_fasta {input.uniprot_fasta} {params.ensembl_fasta_regex} {params.refseq_fasta_regex} --ensembl_fasta {input.ensembl_fasta} --refseq_fasta {input.refseq_fasta} --output {output} &> {log}
        """


rule transcript_to_reference:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        refseq_gtf=rules.prepare_refseq_gtf.output if config["refseq_gtf"] else [],
        ensembl_gtf=rules.prepare_ensembl_gtf.output,
        chr_ref=rules.get_refseq_acc.output if config["refseq_chr_accessions"] else [],
        classification_file=select_sqanti_classification,
        corrected_gtf=rules.run_sqanti.output.gtf, 
        fasta_proteins=rules.run_sqanti.output.fasta_proteins,
        species_db=rules.parse_protein_databases.output
    output:
        protein_assoc=os.path.join("data", prefix, "output", "{db}", "protein_assoc_data.txt"),
        nmd=os.path.join("data", prefix, "output", "{db}", "nmd_data.txt")
    log:
        os.path.join("logs", prefix, "{db}", "transcript_to_reference.log")
    shell:
        """
        scripts/transcript2reference.py --ensembl_gtf {input.ensembl_gtf} \
        --refseq_gtf {input.refseq_gtf} \
        --chr_ref {input.chr_ref} --classification_file {input.classification_file} --corrected_gtf {input.corrected_gtf} \
        --corrected_fasta_proteins {input.fasta_proteins} --output_assoc {output.protein_assoc} --output_nmd {output.nmd} \
        --database {wildcards.db} --species_db {input.species_db} &> {log}    
        """

# LAYERS
rule layer_go:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        classification_file=select_sqanti_classification 
    output:
        _output_layer_db("layer_go", lambda x: config.get("biomart_dataset", []))
    params:
        biomart_host=config.get("biomart_host", []),
        biomart_dataset=config.get("biomart_dataset", []),
    log:
        os.path.join("logs", prefix, "{db}", "layer_go.log")
    shell:
        """
        scripts/layer_go.py --classification_file {input.classification_file} --output {output} --biomart_host {params.biomart_host} --biomart_dataset {params.biomart_dataset} &> {log}
        """


rule layer_interproscan:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        interproscan_file=rules.parse_interproscan.output,
        pfam=rules.prepare_pfam_clan.output,
        classification_file=select_sqanti_classification,
        t=rules.prepare_pfam_clan.input
    output:
       _output_layer_db("layer_interproscan", external_rule=rules.parse_interproscan.output)
    params:
        keep_version=config["transcript_versioned"]
    log:
        os.path.join("logs", prefix, "{db}", "layer_interproscan.log")
    shell:
        """
        (echo {input.t}
        echo {input.pfam}
        scripts/layer_interproscan.py --interproscan_file {input.interproscan_file} --pfam_file {input.pfam} --classification_file {input.classification_file} --keep_version {params.keep_version} --output {output}) &> {log}
        """


rule layer_exons:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        gtf_file=select_gtf,
        chr_ref=rules.get_refseq_acc.output if (config["db"]=="refseq") else []
    output:
        _output_layer_db("layer_exons", external_rule=select_gtf)
    log:
        os.path.join("logs", prefix, "{db}", "layer_exons.log")
    shell:
        """
        scripts/layer_exons.py --gtf_file {input.gtf_file} --chr_ref {input.chr_ref} --output {output} &> {log}
        """


rule layer_junctions:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        junctions_file=lambda x: select_sqanti_output(x).junctions, 
        classification_file=select_sqanti_classification 
    output:
        _output_layer_db("layer_junctions", external_rule=lambda x: select_sqanti_output(x).junctions)
    log:
        os.path.join("logs", prefix, "{db}", "layer_junctions.log")
    shell:
        """
        scripts/layer_junctions.py --junctions_file {input.junctions_file} --classification_file {input.classification_file} --output {output} &> {log}
        """

rule layer_nmd:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        nmd_file=select_nmd_file,
        classification_file=select_sqanti_classification
    output:
        _output_layer_db("layer_nmd", external_rule=select_nmd_file)
    log:
        os.path.join("logs", prefix, "{db}", "layer_nmd.log")
    shell:
        """
        scripts/layer_nmd.py --nmd_file {input.nmd_file} --classification_file {input.classification_file} --output {output} &> {log}
        """


rule layer_reactome:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        classification_file=select_sqanti_classification,
        reactome_file=rules.get_reactome.output
    params:
        biomart_host=config.get("biomart_host", []),
        biomart_dataset=config.get("biomart_dataset", []),
        species=config["species"]
    output:
        _output_layer_db("layer_reactome", external_rule=rules.get_reactome.output)
    log:
        os.path.join("logs", prefix, "{db}", "layer_reactome.log")
    shell:
        """
        scripts/layer_reactome.py --reactome_file {input.reactome_file} --classification_file {input.classification_file} --biomart_host {params.biomart_host} --biomart_dataset {params.biomart_dataset} --species {params.species:q} --output {output} &> {log}
        """


rule layer_repeatmasker:
    conda:
        "../../envs/repeats.yaml"
    input:
        repeatmasker_file=rules.run_repeatmasker.output,
        classification_file=select_sqanti_classification
    output:
        _output_layer_db("layer_repeatmasker", external_rule=rules.run_repeatmasker.output)
    params:
        keep_version=config["transcript_versioned"]
    log:
        os.path.join("logs", prefix, "{db}", "layer_repeatmasker.log")
    shell:
        """
        scripts/layer_repeatmasker.py --keep_version {params.keep_version} --repeatmasker_file {input.repeatmasker_file} --classification_file {input.classification_file} --output {output} &> {log}
        """


rule layer_uniprot:
    input:
        uniprot_file=rules.get_uniprot_phosphosite_annotation.output,
    output:
        _output_layer_db("layer_uniprot", external_rule=rules.get_uniprot_phosphosite_annotation.output)
    log:
        os.path.join("logs", prefix, "{db}", "layer_uniprot.log")
    shell:
        """
        cat {input.uniprot_file} | sort -k1 -k3 -k4 -k5 | uniq | sort -k1 > {output} 2> {log}
        """


rule layer_utrscan:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        utrscan_file=rules.run_utrscan.output,
        classification_file=select_sqanti_classification
    output:
        _output_layer_db("layer_utrscan", external_rule=rules.run_utrscan.output)
    params:
        keep_version=config.get("transcript_versioned", False)
    log:
        os.path.join("logs", prefix, "{db}", "layer_utrscan.log")
    shell:
        """
        scripts/layer_utrscan.py --keep_version {params.keep_version} --utrscan_file {input.utrscan_file} --classification_file {input.classification_file} --output {output} &> {log}
        """


# GET EVERYTHING TOGETHER
rule tappas_annotation:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        transcript_block = [
            rules.layer_utrscan.output,
            rules.layer_repeatmasker.output,
            rules.layer_nmd.output
        ] + config.get("transcript_gtf", []),
        genomic_block = [
            rules.layer_exons.output,
            rules.layer_junctions.output,
        ] + config.get("genomic_gtf", []),
        protein_block = [
            rules.layer_go.output if config["layer_go"] == "si" else [],
            rules.layer_reactome.output if config["reactome"] else [],
            rules.layer_interproscan.output,
            rules.layer_uniprot.output
        ] + config.get("protein_gtf", []),
        classification_file=select_sqanti_classification,
        gene_desc=[],
        protein_assoc=select_prot_assoc
    output:
        os.path.join("data", prefix, f"{species_name}_tappas_{{db}}_annotation_file.gff3")
    log:
        os.path.join("logs", prefix, "{db}", "tappas_annotation.log")
    shell:
        """
        scripts/t2goAnnotationFile.py --classification_file {input.classification_file}  \
         --gene_desc_file {input.gene_desc} --input_transcripts {input.transcript_block} --input_genomic {input.genomic_block} \
         --input_protein  {input.protein_block} --output {output} --gene_desc_file {input.gene_desc} --protein_association {input.protein_assoc} &> {log}
        """


rule renameFeatures:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.tappas_annotation.output
    output:
        os.path.join("data", prefix, f"{species_name}_tappas_{{db}}_annotation_file.gff3_mod") 
    log:
        os.path.join("logs", prefix, "{db}", "renameFeatures.log")
    shell:
        """
        scripts/renameFeatures.py {input} &> {log}
        """

# TODO ASK the user the directory where they want to store the final annotation. Connect with isoannot.sh script
# rule copyToOutDIR:
#     input:
#         rules.renameFeatures.output
#     output:

#     shell:
#         """
#             cp {input} {output}
#         """
