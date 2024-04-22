import os

configfile: "config/generic/config.yaml"
def _output_layer_db(layer_name, external_rule=[], wildcards=None):

    if callable(external_rule):
        external_rule = external_rule({"prefix": config['prefix'], "db": config['db'] })

    if len(config.get(layer_name, [])) or len(external_rule):
        return f"data/{{prefix}}/output/{{db}}/layers/{layer_name}.gtf"
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
        expand("data/{prefix}/{species_name}_tappas_{db}_annotation_file.gff3_mod", prefix = config["prefix"], species_name=config["species_name"], db=config["db"])[0]


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
        expand("data/{{prefix}}/config/refseq/{filename}", filename = os.path.basename(config["refseq_chr_accessions"]))[0]
    params:
        URL=config["refseq_chr_accessions"]
    shell:
        """
        wget -nv -P data/{wildcards.prefix}/config/refseq/ {params.URL}
        """


checkpoint get_refseq_proteins:
    output:
        directory("data/{prefix}/config/refseq/faa/")
    params:
        URL=config["refseq_protein_dir"]
    shell:
        """
        wget -nv -P data/{wildcards.prefix}/config/refseq/faa/ -nd -np -r -nH -l1 -A *protein.faa.gz {params.URL} 
        """


def gather_refseq_proteins(wildcards):
    checkpoint_output = checkpoints.get_refseq_proteins.get(**wildcards).output[0]
    fname_vals = glob_wildcards(os.path.join(checkpoint_output, "{filename}.faa.gz")).filename
    return expand(os.path.join(checkpoint_output, "{filename}.faa.gz"), dir=checkpoint_output, filename=fname_vals)[0]


rule prepare_refseq_proteins:
    input:
       gather_refseq_proteins
    output:
       "data/{prefix}/config/refseq/all_refseq_proteins.faa"
    shell:
       """
       zcat data/{wildcards.prefix}/config/refseq/faa/*.faa.gz > {output}
       """


checkpoint get_refseq_cdna:
    output:
        directory("data/{prefix}/config/refseq/fna/")
    params:
        URL=config["refseq_protein_dir"]
    shell:
        """
        wget -nv -P data/{wildcards.prefix}/config/refseq/fna/ -nd -np -r -nH -l1 -A *rna.fna.gz {params.URL}
        """


def gather_refseq_cdna(wildcards):
    checkpoint_output = checkpoints.get_refseq_cdna.get(**wildcards).output[0]
    fname_vals = glob_wildcards(os.path.join(checkpoint_output, "{filename}.fna.gz")).filename
    return expand(os.path.join(checkpoint_output, "{filename}.fna.gz"), dir=checkpoint_output, filename=fname_vals)[0]


rule prepare_refseq_cdna:
    input:
        gather_refseq_cdna
    output:
        "data/{prefix}/config/refseq/all_refseq_cdna.fna"
    shell:
        """
        zcat data/{wildcards.prefix}/config/refseq/fna/*.fna.gz > {output}
        """

rule get_refseq_gtf:
    output:
        expand("data/{{prefix}}/config/refseq/{filename}", filename = os.path.basename(config["refseq_gtf"]))[0]
    params:
        URL=config["refseq_gtf"]
    shell:
        """
        wget -nv -P data/{wildcards.prefix}/config/refseq/ {params.URL}
        """


rule prepare_refseq_gtf:
    input:
        rules.get_refseq_gtf.output
    output:
        gtf=expand("data/{{prefix}}/config/refseq/{filename}_nopartial_nc.gtf", filename = _remove_extension(config["refseq_gtf"], all=True))[0],
    shell:
        """
        zcat {input}|grep -v -F 'partial=true' | grep NC_0  > {output.gtf}
        """


rule get_ensembl_proteins:
    conda:
        "../../envs/isoannotpy.yaml"   
    output:
        gz=expand("data/{{prefix}}/config/ensembl/{filename}", filename = os.path.basename(config["ensembl_proteins"]))[0],
        fa=expand("data/{{prefix}}/config/ensembl/{filename}", filename = _remove_extension(config["ensembl_proteins"]))[0]

    params:
        URL=config["ensembl_proteins"]
    shell:
        """
        wget -nv -P data/{wildcards.prefix}/config/ensembl/ {params.URL}
        gunzip {output.gz}
        sed -i s/\*//g {output.fa}
        gzip -k {output.fa}
        """


rule get_ensembl_cdna:
    output:
        expand("data/{{prefix}}/config/ensembl/{filename}", filename = os.path.basename(config["ensembl_cdna"]))[0]
    params:
        URL=config["ensembl_cdna"]
    shell:
        """
        wget -nv -P data/{wildcards.prefix}/config/ensembl/ {params.URL}
        """


rule prepare_ensembl_cdna:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.get_ensembl_cdna.output
    output:
        expand("data/{{prefix}}/config/ensembl/{filename}", filename = _remove_extension(config["ensembl_cdna"]))[0]
    shell:
        """
        gunzip -k {input}
        """


rule get_ensembl_reference:
    output:
        expand("data/{{prefix}}/config/ensembl/{filename}", filename = os.path.basename(config["ensembl_reference"]))[0]
    params:
        URL=config["ensembl_reference"]
    shell:
        """
        wget -nv -P data/{wildcards.prefix}/config/ensembl/ {params.URL}
        """

rule prepare_ensembl_reference:
    input:
        rules.get_ensembl_reference.output
    output:
        expand("data/{{prefix}}/config/ensembl/{filename}", filename = _remove_extension(config["ensembl_reference"]))[0]
    shell:
        """
        gunzip -k {input}
        """

rule get_ensembl_gtf:
    output:
        expand("data/{{prefix}}/config/ensembl/{filename}", filename = os.path.basename(config["ensembl_gtf"]))[0]
    params:
        URL=config["ensembl_gtf"]
    shell:
        """
        wget -nv -P data/{wildcards.prefix}/config/ensembl/ {params.URL}
        """


rule prepare_ensembl_gtf:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.get_ensembl_gtf.output
    output:
        expand("data/{{prefix}}/config/ensembl/{filename}", filename = _remove_extension(config["ensembl_gtf"]))[0]
    shell:
        """
        gunzip -k {input}
        """


rule get_pfam_clan:
    output:
        expand("data/global/pfam/{filename}", filename = os.path.basename(config["pfam_clan_url"]))[0]
    params:
        URL=config["pfam_clan_url"]
    shell:
        """
        wget -nv -P data/global/pfam/ {params.URL}
        """


rule prepare_pfam_clan: 
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.get_pfam_clan.output
    output:
        expand("data/global/pfam/{filename}", filename = _remove_extension(config["pfam_clan_url"]))[0]
    shell:
        """
        gunzip -k {input}
        """


rule get_uniprot_data: 
    output:
        expand("data/{{prefix}}/config/uniprot/{uniprot_file}", uniprot_file = [os.path.basename(uniprot_url) for uniprot_url in config["uniprot_dat"] + config["uniprot_fasta"]])
    run:
        for uniprot_url in config["uniprot_dat"] + config["uniprot_fasta"]:
            shell(f"wget -nv -P data/{wildcards.prefix}/config/uniprot/ {uniprot_url}")


rule prepare_uniprot_data:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        [file for file in rules.get_uniprot_data.output if file.endswith(".dat.gz")]
    output:
        "data/{prefix}/config/uniprot/uniprot_parsed.txt"
    shell:
        """
        scripts/uniprot_parse.py --uniprot_files {input} --output {output}
        """


rule get_reactome:
    output:
        expand("data/global/{filename}", filename = os.path.basename(config["reactome"]))[0]
    params:
        URL=config["reactome"]
    shell:
        """
        wget -nv --no-check-certificate -P data/global/ {params.URL}
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
        classification=expand("data/{{prefix}}/config/{{db}}/sqanti_classification.txt", prefix=config["prefix"])[0],
        junctions=expand("data/{{prefix}}/config/{{db}}/sqanti_junctions.txt", prefix=config["prefix"])[0], 
        NMD=expand("data/{{prefix}}/config/{{db}}/sqanti_NMD.txt", prefix=config["prefix"])[0]
    shell:
        """
        scripts/referenceSQANTI.py --gtf_file {input.gtf} --reference_file {input.reference} --chr_ref {input.chr_ref} --database {wildcards.db} --output_classification {output.classification} --output_junctions {output.junctions} --output_nmd {output.NMD}
        """

rule run_gmap_index:
    conda:
        "../../envs/sqanti.yaml"
    input:
        ref_genome=rules.prepare_ensembl_reference.output
    output:
        directory("data/{prefix}/config/ensembl/gmap_index")
    params:
        outdir="data/{prefix}/config/ensembl/",
        index_name="gmap_index"
    shell:
        """
        gmap_build --dir {params.outdir} --genomedb {params.index_name} {input.ref_genome}
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
        corrected_cdna="data/{prefix}/output/{db}/sqanti_corrected.fasta",
        fasta_proteins="data/{prefix}/output/{db}/sqanti_fasta_proteins.fasta",
        gtf="data/{prefix}/output/{db}/sqanti_gtf.gtf",
        classification="data/{prefix}/output/{db}/sqanti_classification.txt",
        junctions="data/{prefix}/output/{db}/sqanti_junctions.txt",
    params:
        outdir="data/{prefix}/output/{db}/",
        out_name="sqanti"
    shell:
        """
        scripts/sqanti/sqanti_qc.py {input.user_fasta} {input.reference_gtf} {input.genome_fasta} -d {params.outdir} -x {input.genome_fasta_index}/gmap_index -o {params.out_name} -n
        """


rule run_utrscan:
    input:
        select_fasta_cdna
    output:
        "data/{prefix}/output/{db}/utrscan.txt"
    shell:
        """
        software/bin/UtrScan -SIGNALLIST -COMMAND=software/bin/UtrSite.Command -INPUT={input} -OUTPUT={output}
        """

rule run_repeatmasker: 
    conda:
        "../../envs/repeats.yaml"
    input:
        select_fasta_cdna
    output:
        expand("data/{{prefix}}/output/{{db}}/repeat_masker/{filename}.out", filename = os.path.basename(select_fasta_cdna(config)))[0]
    params:
        species_name=config["species_name"],
        outdir="data/{prefix}/output/{db}/repeat_masker/"
    shell:
        """
        RepeatMasker {input} -species {params.species_name} -dir {params.outdir}
        """

 
checkpoint run_interproscan:
    conda:
        "../../envs/interpro_java.yaml"
    input:
        select_fasta_proteins
    output:
        directory("data/{prefix}/output/{db}/interproscan/splitProteins/")
    params:
        interproscan_path=config["interproscan_path"]
    shell:
        """
        for i in {input}
        do
        mkdir -p {output} && {params.interproscan_path} -i {input} -d {output} --disable-precalc  -appl Coils,Pfam,MobiDBLite,SignalP_EUK,TMHMM  -f XML -iprlookup
        done
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
        "data/{prefix}/output/{db}/interproscan/interproResultsPfam.tsv"
    shell:
        """
        scripts/parseInterproscanXml.py --interproscan_files {input} --output {output}
        """


rule run_gtftogenepred:
    conda:
        "../../envs/genepred.yaml"
    input:
        select_gtf
    output:
        "data/{prefix}/config/{db}/genePrediction.txt",
    shell:
        """
        gtfToGenePred {input} {output} -genePredExt -allErrors -ignoreGroupsWithoutExons
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
        refseq_gtf = rules.prepare_refseq_gtf.output.gtf if config["refseq_gtf"] else [],
        ensembl_gtf = rules.prepare_ensembl_gtf.output,
        uniprot_parsed = rules.prepare_uniprot_data.output,
        chr_ref = rules.get_refseq_acc.output
    params:
        biomart_host = config.get("biomart_host", [])
    output:
        protein = "data/{prefix}/config/uniprot/uniprot_gcord_proteinGenomic.txt",
        domain = "data/{prefix}/config/uniprot/uniprot_gcord_domainGenomic.txt",
    shell:
        """
        scripts/uniprotPhosphosite_genomicCoordinates.py --uniprot_fasta {input.uniprot_fasta} --refseq_fasta {input.refseq_fasta} --ensembl_fasta {input.ensembl_fasta} --uniprot_parsed {input.uniprot_parsed} --phosphosite_files {input.phosphosite_files} --refseq_gtf {input.refseq_gtf} --ensembl_gtf {input.ensembl_gtf}  --output_protein {output.protein} --output_domain {output.domain} --chr_ref {input.chr_ref} --biomart_host {params.biomart_host}
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
        expand("data/{{prefix}}/output/{{db}}/uniprot_Phosphosite_info.txt")[0]
    shell:           
        """
        scripts/uniprotPhosphosite_annotation.py --orf_fasta {input.fasta_orf} --classification_file {input.classification} --genepred_file {input.gene_prediction} --uniprotmotif_file {input.motif_info} --chr_ref {input.chr_ref} --protein_association {input.protein_assoc} --keep_version {params.keep_version} --db {wildcards.db} --biomart_host {params.biomart_host} --output {output}
        """


rule parse_protein_databases:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        uniprot_fasta = [file for file in rules.get_uniprot_data.output if file.endswith(".fasta.gz")],
        refseq_fasta=rules.prepare_refseq_proteins.output,
        ensembl_fasta=rules.get_ensembl_proteins.output.fa,
    output:
        "data/{prefix}/config/parsed_databases.json"
    params:
        ensembl_fasta_regex=_optional_param("ensembl_fasta_regex", "--ensembl_fasta_regex"),
        refseq_fasta_regex=_optional_param("refseq_fasta_regex", "--refseq_fasta_regex"),
    shell:
        """
        scripts/parse_protein_databases.py --uniprot_fasta {input.uniprot_fasta} {params.ensembl_fasta_regex} {params.refseq_fasta_regex} --ensembl_fasta {input.ensembl_fasta} --refseq_fasta {input.refseq_fasta} --output {output}
        """


rule transcript_to_reference:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        refseq_gtf=rules.prepare_refseq_gtf.output.gtf if config["refseq_gtf"] else [],
        ensembl_gtf=rules.prepare_ensembl_gtf.output,
        chr_ref=rules.get_refseq_acc.output if config["refseq_chr_accessions"] else [],
        classification_file=select_sqanti_classification,
        corrected_gtf=rules.run_sqanti.output.gtf, 
        fasta_proteins=rules.run_sqanti.output.fasta_proteins,
        species_db=rules.parse_protein_databases.output
    output:
        protein_assoc="data/{prefix}/output/{db}/protein_assoc_data.txt",
        nmd="data/{prefix}/output/{db}/nmd_data.txt"
    shell:
        """
        scripts/transcript2reference.py --ensembl_gtf {input.ensembl_gtf} \
        --refseq_gtf {input.refseq_gtf} \
        --chr_ref {input.chr_ref} --classification_file {input.classification_file} --corrected_gtf {input.corrected_gtf} \
        --corrected_fasta_proteins {input.fasta_proteins} --output_assoc {output.protein_assoc} --output_nmd {output.nmd} \
        --database {wildcards.db} --species_db {input.species_db}
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
    shell:
        """
        scripts/layer_go.py --classification_file {input.classification_file} --output {output} --biomart_host {params.biomart_host} --biomart_dataset {params.biomart_dataset} """


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
    shell:
        """
        echo {input.t}
        echo {input.pfam}
        scripts/layer_interproscan.py --interproscan_file {input.interproscan_file} --pfam_file {input.pfam} --classification_file {input.classification_file} --keep_version {params.keep_version} --output {output}
        """


rule layer_exons:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        gtf_file=select_gtf,
        chr_ref=rules.get_refseq_acc.output if (config["db"]=="refseq") else []
    output:
        _output_layer_db("layer_exons", external_rule=select_gtf)
    shell:
        """
        scripts/layer_exons.py --gtf_file {input.gtf_file} --chr_ref {input.chr_ref} --output {output}
        """


rule layer_junctions:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        junctions_file=lambda x: select_sqanti_output(x).junctions, 
        classification_file=select_sqanti_classification 
    output:
        _output_layer_db("layer_junctions", external_rule=lambda x: select_sqanti_output(x).junctions)
    shell:
        """
        scripts/layer_junctions.py --junctions_file {input.junctions_file} --classification_file {input.classification_file} --output {output}
        """

rule layer_nmd:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        nmd_file=select_nmd_file,
        classification_file=select_sqanti_classification
    output:
        _output_layer_db("layer_nmd", external_rule=select_nmd_file)
    shell:
        """
        scripts/layer_nmd.py --nmd_file {input.nmd_file} --classification_file {input.classification_file} --output {output}
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
    shell:
        """
        scripts/layer_reactome.py --reactome_file {input.reactome_file} --classification_file {input.classification_file} --biomart_host {params.biomart_host} --biomart_dataset {params.biomart_dataset} --species {params.species:q} --output {output}
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
    shell:
        """
        scripts/layer_repeatmasker.py --keep_version {params.keep_version} --repeatmasker_file {input.repeatmasker_file} --classification_file {input.classification_file} --output {output}
        """


rule layer_uniprot:
    input:
        uniprot_file=rules.get_uniprot_phosphosite_annotation.output,
    output:
        _output_layer_db("layer_uniprot", external_rule=rules.get_uniprot_phosphosite_annotation.output)
    shell:
        """
        cat {input.uniprot_file} | sort -k1 -k3 -k4 -k5 | uniq | sort -k1 > {output}
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
    shell:
        """
        scripts/layer_utrscan.py --keep_version {params.keep_version} --utrscan_file {input.utrscan_file} --classification_file {input.classification_file} --output {output}
        """


# GET EVERYTHING TOGETHER
rule tappas_annotation:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        transcript_block = [
            # rules.layer_utrscan.output,
            # rules.layer_repeatmasker.output,   # "data/{prefix}/output/{db}/layers/layer_repeatmasker.gtf",
            # rules.layer_nmd.output
        ] + config.get("transcript_gtf", []),
        genomic_block = [
            # rules.layer_exons.output,
            # rules.layer_junctions.output,
        ] + config.get("genomic_gtf", []),
        protein_block = [
            rules.layer_go.output if config["layer_go"] == "si" else [],
            # rules.layer_reactome.output if config["reactome"] else [],
            # rules.layer_interproscan.output,
            # "data/{prefix}/output/{db}/layers/layer_uniprot.gtf"
        ] + config.get("protein_gtf", []),
        classification_file=select_sqanti_classification,
        gene_desc=[],
        protein_assoc=select_prot_assoc
    output:
        expand("data/{{prefix}}/{species_name}_tappas_{{db}}_annotation_file.gff3", species_name = config["species_name"])[0] 
    shell:
        """
        scripts/t2goAnnotationFile.py --classification_file {input.classification_file}  \
         --gene_desc_file {input.gene_desc} --input_transcripts {input.transcript_block} --input_genomic {input.genomic_block} \
         --input_protein  {input.protein_block} --output {output} --gene_desc_file {input.gene_desc} --protein_association {input.protein_assoc}
        """


rule renameFeatures:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.tappas_annotation.output
    output:
        expand("data/{{prefix}}/{species}_tappas_{{db}}_annotation_file.gff3_mod", species = config["species_name"])[0] 
    shell:
        """
        scripts/renameFeatures.py {input}
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
