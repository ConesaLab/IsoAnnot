import os

configfile: "config/generic/config.yaml"

def is_fasta_input():
    db = config.get("db", "ensembl")
    if db in ["ensembl", "refseq"]:
        return True
    return config.get("fasta_cdna", None) is not None

def select_gmap_index(wildcards):
    if not is_fasta_input():
        return []
    return rules.run_gmap_index.output

prefix = config["prefix"]
db = config["db"]
species_name = config["species_name"]
path_output = config["path_output"]
nls_model = config.get("nls_model", None)
nls_chunks = 50


def is_feature_enabled(key, default=True):
    val = config.get(key, None)
    if val is None:
        return default
    if isinstance(val, bool):
        return val
    if isinstance(val, str):
        return val.lower() not in ["false", "no", "0", "off"]
    return bool(val)


def _output_layer_db(layer_name, external_rule=[], wildcards=None):
    if layer_name == "layer_repeatmasker" and not is_feature_enabled("repeat_masker", is_feature_enabled("repeatmasker", True)):
        return []
    if callable(external_rule):
        external_rule = external_rule({"prefix": prefix, "db": db })
    if len(config.get(layer_name, [])) or len(external_rule):
        return os.path.join(path_output, "data", prefix, "output", "{db}", "layers", f"{layer_name}.gtf")
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

def get_sqanti_extra_params(wildcards):
    if config.get("fasta_cdna"):
        return "--fasta"
    return ""

rule all:
    input:
        expand("{path_output}/data/{prefix}/{species_name}_tappas_{db}_annotation_file.gff3_mod",
               path_output=path_output,
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
        os.path.join(path_output, "data",prefix,"config","refseq", os.path.basename(config["refseq_chr_accessions"])) 
    params: 
        URL=config["refseq_chr_accessions"] 
    log: 
        os.path.join(path_output, "logs", prefix, "get_refseq_acc.log")
    shell: 
        """ 
        wget -nv -P {path_output}/data/{prefix}/config/refseq/ {params.URL} &> {log}
        """

rule get_refseq_proteins:
    output:
        os.path.join(path_output, "data",prefix,"config","refseq", os.path.basename(config["refseq_proteins"]))
    params:
        URL=config["refseq_proteins"]
    log:
        os.path.join(path_output, "logs", prefix, "get_refseq_proteins.log")
    shell:
        """
        wget -nv -P {path_output}/data/{prefix}/config/refseq/ {params.URL} &> {log}
        """

rule prepare_refseq_proteins:
    input:
        rules.get_refseq_proteins.output
    output:
        os.path.join(path_output, "data", prefix, "config", "refseq", f"{_remove_extension(config['refseq_proteins'])}")
    log:
        os.path.join(path_output, "logs", prefix, "prepare_refseq_proteins.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """

rule get_refseq_cdna:
    output:
        os.path.join(path_output, "data",prefix,"config","refseq", os.path.basename(config["refseq_cdna"]))
    params:
        URL=config["refseq_cdna"]
    log:
        os.path.join(path_output, "logs", prefix, "get_refseq_cdna.log")
    shell:
        """
        wget -nv -P {path_output}/data/{prefix}/config/refseq/ {params.URL} &> {log}
        """

rule prepare_refseq_cdna:
    input:
        rules.get_refseq_cdna.output
    output:
        fa=os.path.join(path_output, "data", prefix, "config", "refseq", f"{_remove_extension(config['refseq_cdna'])}")
    log:
        os.path.join(path_output, "logs", prefix, "prepare_refseq_cdna.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """

rule get_refseq_gtf:
    output:
        os.path.join(path_output, "data", prefix, "config", "refseq", os.path.basename(config["refseq_gtf"]))
    params:
        URL=config["refseq_gtf"]
    log:
        os.path.join(path_output, "logs", prefix, "get_refseq_gtf.log")
    shell:
        """
        wget -nv -P {path_output}/data/{prefix}/config/refseq/ {params.URL} &> {log}
        """


rule prepare_refseq_gtf:
    input:
        rules.get_refseq_gtf.output
    output:
        os.path.join(path_output, "data", prefix, "config", "refseq", f"{_remove_extension(config['refseq_gtf'], all=True)}_nopartial_nc.gtf")
    log:
        os.path.join(path_output, "logs", prefix, "prepare_refseq_gtf.log")
    shell:
        """
        zcat {input}|grep -v -F 'partial=true' | grep '^NC_' > {output} 2> {log}
        """


rule get_ensembl_proteins:
    output:
        gz=os.path.join(path_output, "data", prefix, "config", "ensembl", os.path.basename(config["ensembl_proteins"])),
        fa=os.path.join(path_output, "data", prefix, "config", "ensembl", _remove_extension(config["ensembl_proteins"]))
    params:
        URL=config["ensembl_proteins"]
    log:
        os.path.join(path_output, "logs", prefix, "get_ensembl_proteins.log")
    shell:
        """
        (wget -nv -P {path_output}/data/{prefix}/config/ensembl/ {params.URL}
        gunzip -c {output.gz} > {output.fa}.tmp
        sed s/\\*//g {output.fa}.tmp > {output.fa}
        rm {output.fa}.tmp) 2> {log}
        """
     

rule get_ensembl_cdna:
    output:
        os.path.join(path_output, "data", prefix, "config", "ensembl", os.path.basename(config["ensembl_cdna"]))
    params:
        URL=config["ensembl_cdna"]
    log:
        os.path.join(path_output, "logs", prefix, "get_ensembl_cdna.log")
    shell:
        """
        wget -nv -P {path_output}/data/{prefix}/config/ensembl/ {params.URL}
	"""


rule prepare_ensembl_cdna:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.get_ensembl_cdna.output
    output:
        fa=os.path.join(path_output, "data", prefix, "config", "ensembl", _remove_extension(config["ensembl_cdna"]))
    log:
        os.path.join(path_output, "logs", prefix, "prepare_ensembl_cdna.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """


rule get_ensembl_reference:
    output:
        os.path.join(path_output, "data", prefix, "config", "ensembl", os.path.basename(config["ensembl_reference"]))
    params:
        URL=config["ensembl_reference"]
    log:
        os.path.join(path_output, "logs", prefix, "get_ensembl_reference.log")
    shell:
        """
        wget -nv -P {path_output}/data/{prefix}/config/ensembl/ {params.URL}
        """

rule prepare_ensembl_reference:
    input:
        rules.get_ensembl_reference.output
    output:
        os.path.join(path_output, "data", prefix, "config", "ensembl", _remove_extension(config["ensembl_reference"]))
    log:
        os.path.join(path_output, "logs", prefix, "prepare_ensembl_reference.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """

rule get_ensembl_gtf:
    output:
        os.path.join(path_output, "data", prefix, "config", "ensembl", os.path.basename(config["ensembl_gtf"]))
    params:
        URL=config["ensembl_gtf"]
    log:
        os.path.join(path_output, "logs", prefix, "get_ensembl_gtf.log")
    shell:
        """
        wget -nv -P {path_output}/data/{prefix}/config/ensembl/ {params.URL}
        """


rule prepare_ensembl_gtf:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.get_ensembl_gtf.output
    output:
        os.path.join(path_output, "data", prefix, "config", "ensembl", _remove_extension(config["ensembl_gtf"]))
    log:
        os.path.join(path_output, "logs", prefix, "prepare_ensembl_gtf.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """

rule check_chromosome_consistency:
    input:
        user_input = select_user_cdna,
        ref_gtf = select_reference_gtf,
        ref_fasta = rules.prepare_ensembl_reference.output
    output:
        touch(os.path.join(path_output, "data", prefix, "config", "check_chromosomes.done"))
    shell:
        """
        python3 scripts/check_chromosomes.py \
            --user_input {input.user_input} \
            --ref_gtf {input.ref_gtf} \
            --ref_fasta {input.ref_fasta}
        """

rule get_pfam_clan:
    output:
        os.path.join(path_output, "data", "global", "pfam", os.path.basename(config["pfam_clan_url"]))
    params:
        URL=config["pfam_clan_url"]
    log:
        os.path.join(path_output, "logs", prefix, "get_pfam_clan.log")
    shell:
        """
        wget -nv -P {path_output}/data/global/pfam/ {params.URL} &> {log}
        """

rule prepare_pfam_clan: 
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.get_pfam_clan.output
    output:
        os.path.join(path_output, "data", "global", "pfam", _remove_extension(config["pfam_clan_url"]))
    log:
        os.path.join(path_output, "logs", prefix, "prepare_pfam_clan.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """


rule get_uniprot_data: 
    resources:
        n_downloads=1
    output:
        [os.path.join(path_output, "data",prefix,"config","uniprot",os.path.basename(uniprot_url)) for uniprot_url in config["uniprot_dat"] + config["uniprot_fasta"]]
    log:
        [os.path.join(path_output, "logs", prefix, "get_uniprot_data", f"{_remove_extension(uniprot_url)}.log") for uniprot_url in config["uniprot_dat"] + config["uniprot_fasta"]]
    run:
        for uniprot_url, out_file, log_file in zip(config["uniprot_dat"] + config["uniprot_fasta"], output, log):
            shell(f"wget -nv -P {path_output}/data/{prefix}/config/uniprot/ {uniprot_url} &> {log_file}")


rule prepare_uniprot_data:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        [file for file in rules.get_uniprot_data.output if file.endswith(".dat.gz")]
    output:
        os.path.join(path_output, "data", prefix, "config", "uniprot", "uniprot_parsed.txt")
    log:
        os.path.join(path_output, "logs", prefix, "prepare_uniprot_data.log")
    shell:
        """
        scripts/uniprot_parse.py --uniprot_files {input} --output {output} &> {log}
        """


rule get_reactome:
    conda:
        "../../envs/isoannotpy.yaml" 
    output:
        os.path.join(path_output, "data", "global", os.path.basename(config["reactome"]))
    params:
        URL=config["reactome"]
    log:
        os.path.join(path_output, "logs", prefix, "get_reactome.log")
    shell:
        """
        curl -L -k -o {output} {params.URL} &> {log}
        """

rule get_mirwalk:
    output:
        os.path.join(path_output, "data", prefix, "config", "mirna", "{region}.zip")
    log:
        os.path.join(path_output, "logs", prefix, "get_mirwalk_{region}.log")
    run:
        target_url = config["mirwalk_urls"].get(wildcards.region)

        if not target_url:
            error_msg = (
                f"FATAL ERROR: You are attempting to analyze region '{wildcards.region}', "
                f"but a valid URL is not defined in config['mirwalk_urls']."
            )
            raise ValueError(error_msg)

        shell("wget -O {output} {target_url} &> {log}")


rule prepare_mirwalk:
    input:
        files = expand(
            os.path.join(path_output, "data", prefix, "config", "mirna", "{region}.zip"),
            region=config.get("mirna_regions_to_use", [])
        )
    output:
        merged = os.path.join(path_output, "data", prefix, "config", "mirna", "mirwalk_merged.txt")
    log:
        os.path.join(path_output, "logs", prefix, "prepare_mirwalk.log")
    shell:
        """
        > {output.merged}
        if [ -z "{input.files}" ]; then
            echo "WARNING: No regions defined. Output empty." > {log}
        else
            for file in {input.files}; do
                echo "Processing $file..." >> {log}
                unzip -p "$file" >> {output.merged} 2>> {log}
            done
        fi
        """
	

rule get_rna_fasta_mirwalk:
    output:
        os.path.join(path_output, "data", prefix, "config", "mirna", "mirwalk_rna_reference.fna.gz")
    params:
        URL = config.get("rna_fasta_mirwalk", "")
    log:
        os.path.join(path_output, "logs", prefix, "get_rna_fasta_mirwalk.log")
    shell:
        """
        wget -nv -O {output} {params.URL} &> {log}
        """

rule prepare_rna_fasta_mirwalk:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.get_rna_fasta_mirwalk.output
    output:
         os.path.join(path_output, "data", prefix, "config", "mirna", "mirwalk_rna_reference.fna")
    log:
         os.path.join(path_output, "logs", prefix, "prepare_rna_fasta_mirwalk.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """

rule get_gtf_mirwalk:
    output:
        os.path.join(path_output, "data", prefix, "config", "mirna", "mirwalk_genomic_reference.gtf.gz")
    params:
        URL = config.get("gtf_mirwalk", "")
    log:
        os.path.join(path_output, "logs", prefix, "get_gtf_mirwalk.log")
    shell:
        """
        wget -nv -O {output} {params.URL} &> {log}
        """

rule prepare_gtf_mirwalk:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.get_gtf_mirwalk.output
    output:
        os.path.join(path_output, "data", prefix, "config", "mirna", "mirwalk_genomic_reference.gtf")
    log:
        os.path.join(path_output, "logs", prefix, "prepare_gtf_mirwalk.log")
    shell:
        """
        gunzip -k {input} &> {log}
        """

# RUN

rule run_refsqanti:
    conda:
        "../../envs/sqanti3.yaml"
    input:
        gtf=select_reference_gtf,
        reference=rules.prepare_ensembl_reference.output,
        chr_ref=rules.get_refseq_acc.output  # we use this even in ensembl mode to filter only chromosome sequences (avoid MT)
    output:
        classification=os.path.join(path_output, "data", prefix, "config", "{db}", "sqanti_classification.txt"),
        junctions=os.path.join(path_output, "data", prefix, "config", "{db}", "sqanti_junctions.txt"),
        NMD=os.path.join(path_output, "data", prefix, "config", "{db}", "sqanti_NMD.txt")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "run_refsqanti.log")
    shell:
        """
        scripts/referenceSQANTI.py --gtf_file {input.gtf} --reference_file {input.reference} --chr_ref {input.chr_ref} --database {wildcards.db} --output_classification {output.classification} --output_junctions {output.junctions} --output_nmd {output.NMD} &> {log}
        """

rule run_gmap_index:
    conda:
        "../../envs/sqanti3.yaml"
    input:
        ref_genome=rules.prepare_ensembl_reference.output
    output:
        directory(os.path.join(path_output, "data", prefix, "config", "ensembl", "gmap_index"))
    params:
        outdir=os.path.join(path_output, "data",prefix,"config","ensembl"),
        index_name="gmap_index"
    log:
        os.path.join(path_output, "logs", prefix, "run_gmap_index.log")
    shell:
        """
        gmap_build -D {params.outdir} -d {params.index_name} {input.ref_genome} &> {log}
        """
rule run_sqanti:
    conda:
        "../../envs/sqanti3.yaml"
    input:
        user_cdna=select_cdna_with_warnings if "select_cdna_with_warnings" in globals() else select_user_cdna,
        reference_gtf=select_reference_gtf,
        genome_fasta=rules.prepare_ensembl_reference.output,
        genome_fasta_index=select_gmap_index,
        check = rules.check_chromosome_consistency.output
    output:
        corrected_cdna=os.path.join(path_output, "data", prefix, "output", "{db}", "sqanti_corrected.fasta"),
        fasta_proteins=os.path.join(path_output, "data", prefix, "output", "{db}", "sqanti_corrected.faa"),
        gtf=os.path.join(path_output, "data", prefix, "output", "{db}", "sqanti_corrected.gtf"),
        classification=os.path.join(path_output, "data", prefix, "output", "{db}", "sqanti_classification.txt"),
        junctions=os.path.join(path_output, "data", prefix, "output", "{db}", "sqanti_junctions.txt"),
    params:
        outdir=os.path.join(path_output, "data",prefix,"output","{db}"),
        out_name="sqanti",
        extra_flag = get_sqanti_extra_params,
        gmap_option=lambda wildcards: f"-x {os.path.join(path_output, 'data', prefix, 'config', 'ensembl', 'gmap_index', 'gmap_index')}" if is_fasta_input() else ""
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "run_sqanti.log")
    shell:
        """
        sqanti3_qc.py --isoforms {input.user_cdna} --refGTF {input.reference_gtf} --refFasta {input.genome_fasta} \
            -d {params.outdir} -o {params.out_name} \
            --include_ORF \
            {params.gmap_option} {params.extra_flag} &> {log}
        """

rule clean_sqanti_proteins:
    input:
        rules.run_sqanti.output.fasta_proteins
    output:
        os.path.join(path_output, "data", prefix, "output", "{db}", "final_sqanti_corrected.faa")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "clean_sqanti_proteins.log")
    shell:
        """
        sed s/\\*//g {input} > {output} 2> {log}
        """

rule run_utrscan:
    input:
        select_fasta_cdna
    output:
        os.path.join(path_output, "data", prefix, "output", "{db}", "utrscan.txt")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "run_utrscan.log")
    params:
        utrscan_bin = os.path.abspath("software/bin/UtrScan"),
        utrsite_cmd = os.path.abspath("software/bin/UtrSite.Command")
    shell:
        """
        {params.utrscan_bin} -SIGNALLIST -COMMAND={params.utrsite_cmd} -INPUT={input} -OUTPUT={output} &> {log}
        """

rule run_repeatmasker:
    conda:
        "../../envs/repeats.yaml"
    input:
        select_fasta_cdna
    output:
        os.path.join(path_output, "data", prefix, "output", "{db}", "repeat_masker", f"{os.path.basename(select_fasta_cdna(config))}.out")
    params:
        species=config["species"],
        outdir=os.path.join(path_output, "data", prefix, "output", "{db}", "repeat_masker/"),
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "run_repeatmasker.log")
    shell:
        r"""
        mkdir -p {params.outdir}
        cd {params.outdir}

        # Prepare libraries
        LIBDIR="$CONDA_PREFIX/lib"
        if [ ! -e "$LIBDIR/libnsl.so.1" ]; then
            ln -sf "$LIBDIR/libnsl.so.3" "$LIBDIR/libnsl.so.1"  #A symbolic link is created; there is no version that installs the required dependency
        fi
        export LD_LIBRARY_PATH="$LIBDIR:$LD_LIBRARY_PATH"

        # Create a mapping file and a FASTA with short IDs (MD5 Hashes)
        # We only apply hashing if the ID length is > 50 characters to avoid RepeatMasker errors
        tmp_fasta="tmp_hashed.fasta"
        mapping_file="id_mapping.tsv"

        python3 -c "
import sys, hashlib
with open('{input}', 'r') as f, open('$tmp_fasta', 'w') as out, open('$mapping_file', 'w') as m:
    for line in f:
        if line.startswith('>'):
            original_id = line[1:].strip().split()[0]
            if len(original_id) > 50:
                new_id = hashlib.md5(original_id.encode()).hexdigest()
                m.write(f'{{new_id}}\\t{{original_id}}\\n')
                out.write(f'>{{new_id}}\\n')
            else:
                out.write(line)
        else:
            out.write(line)
"

        # Run RepeatMasker
        RepeatMasker $tmp_fasta -species "{params.species}" -dir . &> {log}

        # Restore original IDs in the .out file
        # Use a small Python script to replace the temporary hashes with the original long names
        python3 -c "
import os
if os.path.exists(tmp_fasta):
mapping = dict(line.strip().split('\\t') for line in open('$mapping_file'))
with open('$tmp_fasta.out', 'r') as f_in, open('final_corrected.out', 'w') as f_out:
    for line in f_in:
        for short_id, long_id in mapping.items():
            line = line.replace(short_id, long_id)
        f_out.write(line)
"
        # 5. Cleanup
        mv final_corrected.out {output}
        """

rule filter_interactions:
    input:
        mirwalk=rules.prepare_mirwalk.output
    params:
        species_name=config["species"],
        mirbase="data/global/miRNA/miRNA.dat",
        score=config.get("mirna_db_evidence_score_threshold")
    output:
        os.path.join(path_output, "data", prefix, "config", "mirna", "filter_interactions.txt")
    log:
        os.path.join(path_output, "logs", prefix, "filter_interactions.log")
    shell:
        """
        scripts/filter_mirna_bs.py --mirbase_file {params.mirbase} --mirwalk_file {input.mirwalk} --species {params.species_name:q} --score {params.score} --mirwalk_output {output} &> {log}
        """
 
rule run_mirwalk2gen:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        mirwalk=rules.filter_interactions.output,
        fasta=rules.prepare_rna_fasta_mirwalk.output,
        gtf=rules.prepare_gtf_mirwalk.output,
        chr_ref=rules.get_refseq_acc.output
    output:
        os.path.join(path_output, "data", prefix, "output", "{db}", "interactions_gc.txt")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "run_mirwalk2gen.log")
    shell:
        """
        scripts/mirna_bs_genomic_coord.py --mirwalk_file {input.mirwalk} --refseq_fasta {input.fasta} --refseq_gtf {input.gtf} --mirna_output {output} --chr_ref {input.chr_ref} &> {log}
        """


checkpoint run_interproscan:
    conda:
        "../../envs/interpro_java.yaml"
    input:
        select_fasta_proteins
    output:
        directory(os.path.join(path_output, "data", prefix, "output", "{db}", "interproscan", "splitProteins"))  
    params:
        interproscan_path=config["interproscan_path"]
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "run_interproscan.log")
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
        "../../envs/sqanti3.yaml"
    input:
        gather_interproscan
    output:
        os.path.join(path_output, "data", prefix, "output", "{db}", "interproscan", "interproResultsPfam.tsv")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "parse_interproscan.log")
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
        os.path.join(path_output, "data", prefix, "config", "{db}", "genePrediction.txt")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "run_gtftogenepred.log")
    shell:
        """
        gtfToGenePred {input} {output} -genePredExt -allErrors -ignoreGroupsWithoutExons &> {log}
        """

rule get_mirna_bs_annotation:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        fasta=select_fasta_cdna,
        gene_prediction=rules.run_gtftogenepred.output,
        mirna_bs=rules.run_mirwalk2gen.output,
        chr_ref=rules.get_refseq_acc.output
    params:
        db=config.get("db")
    output:
        os.path.join(path_output, "data", prefix, "output", "{db}", "mirna_bs_annotation.txt")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "get_mirna_bs_annotation.log")
    shell:
        """
        scripts/get_mirna_bs_annotation.py --chr_ref {input.chr_ref} --genepred {input.gene_prediction} --isoform_fasta {input.fasta} --mirwalk_genomic {input.mirna_bs} --output {output} --db {params.db} &> {log}
        """

rule get_genomic_coordinates:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        uniprot_fasta = [file for file in rules.get_uniprot_data.output if file.endswith(".fasta.gz")],
        refseq_fasta = rules.prepare_refseq_proteins.output,
        ensembl_fasta = rules.get_ensembl_proteins.output.fa,
        phosphosite_files = [
            "data/global/PSP_data/Acetylation_site_dataset",
            "data/global/PSP_data/Methylation_site_dataset",
            "data/global/PSP_data/O-GalNAc_site_dataset",
            "data/global/PSP_data/O-GlcNAc_site_dataset",
            "data/global/PSP_data/Phosphorylation_site_dataset",
            "data/global/PSP_data/Sumoylation_site_dataset",
            "data/global/PSP_data/Ubiquitination_site_dataset"],
        refseq_gtf = rules.prepare_refseq_gtf.output if config.get("refseq_gtf") else [],
        ensembl_gtf = rules.prepare_ensembl_gtf.output,
        uniprot_parsed = rules.prepare_uniprot_data.output,
        chr_ref = rules.get_refseq_acc.output
    output:
        protein=os.path.join(path_output, "data", prefix, "config", "uniprot", "uniprot_gcord_proteinGenomic.txt"),
        domain=os.path.join(path_output, "data", prefix, "config", "uniprot", "uniprot_gcord_domainGenomic.txt"),
    params:
        biomart_host = config.get("biomart_host", [])
    log:
        os.path.join(path_output, "logs", prefix, "get_genomic_coordinates.log")
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
        os.path.join(path_output, "data", prefix, "output", "{db}", "uniprot_Phosphosite_info.txt")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "get_uniprot_phosphosite_annotation.log")
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
        os.path.join(path_output, "data", prefix, "config", "parsed_databases.json")
    params:
        ensembl_fasta_regex=_optional_param("ensembl_fasta_regex", "--ensembl_fasta_regex"),
        refseq_fasta_regex=_optional_param("refseq_fasta_regex", "--refseq_fasta_regex"),
    log:
        os.path.join(path_output, "logs", prefix, "parse_protein_databases.log")
    shell:
        """
        scripts/parse_protein_databases.py --uniprot_fasta {input.uniprot_fasta} {params.ensembl_fasta_regex} {params.refseq_fasta_regex} --ensembl_fasta {input.ensembl_fasta} --refseq_fasta {input.refseq_fasta} --output {output} &> {log}
        """


rule transcript_to_reference:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        refseq_gtf=rules.prepare_refseq_gtf.output if config.get("refseq_gtf") else [],
        ensembl_gtf=rules.prepare_ensembl_gtf.output,
        chr_ref=rules.get_refseq_acc.output if config.get("refseq_chr_accessions") else [],
        classification_file=select_sqanti_classification,
        corrected_gtf=rules.run_sqanti.output.gtf, 
        fasta_proteins=rules.run_sqanti.output.fasta_proteins,
        species_db=rules.parse_protein_databases.output
    output:
        protein_assoc=os.path.join(path_output, "data", prefix, "output", "{db}", "protein_assoc_data.txt"),
        nmd=os.path.join(path_output, "data", prefix, "output", "{db}", "nmd_data.txt")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "transcript_to_reference.log")
    shell:
        """
        scripts/transcript2reference.py --ensembl_gtf {input.ensembl_gtf} \
        --refseq_gtf {input.refseq_gtf} \
        --chr_ref {input.chr_ref} --classification_file {input.classification_file} --corrected_gtf {input.corrected_gtf} \
        --corrected_fasta_proteins {input.fasta_proteins} --output_assoc {output.protein_assoc} --output_nmd {output.nmd} \
        --database {wildcards.db} --species_db {input.species_db} &> {log}    
        """

rule filter_nls:
    conda:
        "../../envs/isoannotpy.yaml"
    input: 
        select_fasta_proteins
    output: 
        os.path.join(path_output,"data",prefix,"output","{db}","nls","nls_filtered_proteins.fa")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "filter_nls.log")
    shell: 
        """
        scripts/nls_filter.py --input {input} --output {output} &> {log}
        """

rule nls_deduplicate:
    conda: 
        "../../envs/isoannotpy.yaml"
    input: 
        rules.filter_nls.output
    output:
        fasta = os.path.join(path_output,"data",prefix,"output","{db}","nls","unique_proteins.fa"),
        mapping = os.path.join(path_output,"data",prefix,"output","{db}","nls","protein_mapping.tsv")
    log: 
        os.path.join(path_output, "logs", prefix, "{db}", "nls_deduplicate.log")
    shell:
        """
        scripts/deduplicate_proteins.py --input {input} \
            --output_fasta {output.fasta} --output_mapping {output.mapping} &> {log}
        """

rule split_proteins:
    input: 
        rules.nls_deduplicate.output.fasta
    output: 
        expand(os.path.join(path_output,"data",prefix,"output","{{db}}","nls", "chunks_temp", "chunk_{n}.fa"), n=range(nls_chunks))
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "split_proteins.log")
    shell:
        """
        mkdir -p $(dirname {output[0]})
        awk 'BEGIN {{RS=">"; FS="\\n"}} \
             NR>1 {{ \
                out_file = "{path_output}/data/{prefix}/output/{wildcards.db}/nls/chunks_temp/chunk_" (i++ % {nls_chunks}) ".fa"; \
                print ">"$0 > out_file; \
             }}' {input} 2> {log}
        """

rule run_nucimport:
    input: 
        os.path.join(path_output, "data", prefix, "output", "{db}", "nls", "chunks_temp", "chunk_{n}.fa")
    output: 
        os.path.join(path_output,"data",prefix,"output","{db}","nls","tmp", "output_chunk_{n}.txt")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "run_nucimport", "run_nucimport_chunk_{n}.log")
    params:
        jar_dir = "software/NucImport",
        jar_name = "NucImportMay2012.jar"
    shell:
        """
        cd {params.jar_dir}
        java -jar {params.jar_name} {input} {nls_model} Mouse ID=F > {output} 2> {log}
        """

rule merge_nls_chunks:
    input: 
        expand(os.path.join(path_output,"data",prefix,"output","{{db}}","nls","tmp", "output_chunk_{n}.txt"), n=range(nls_chunks))
    output: 
        os.path.join(path_output,"data",prefix,"output","{db}","nls","output_merged.txt")
    log:
        os.path.join(path_output,"logs",prefix,"{db}","nls_chunks_merged.log")
    shell:
        """
        awk "NR == FNR || (FNR > 3 && !/^Protein/ && !/^\\*/)" {input} > {output} 2> {log}
        """

rule parse_nls:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.merge_nls_chunks.output
    output:
        os.path.join(path_output,"data",prefix,"output","{db}","nls","nls_parsed.tsv")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "parse_nls.log")
    params:
        t_imp = config.get("nls_threshold_import", 0.7),
        t_cnls = config.get("nls_threshold_cnls", 0.3)
    shell:
        """
        scripts/parse_nls.py --input {input} \
            --threshold_imp {params.t_imp} --threshold_cnls {params.t_cnls} \
            --output {output} &> {log}
        """

rule expand_nls:
    conda: 
        "../../envs/isoannotpy.yaml"
    input:
        parsed_tsv = rules.parse_nls.output,
        mapping = rules.nls_deduplicate.output.mapping
    output:
        os.path.join(path_output,"data",prefix,"output","{db}","nls","nls_final_expanded.tsv")
    log:
       os.path.join(path_output, "logs", prefix, "{db}", "expand_nls.log")
    shell:
        """
        scripts/expand_nls_annotation.py --parsed {input.parsed_tsv} \
            --mapping {input.mapping} --output {output}
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
        db=config.get("db")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "layer_go.log")
    shell:
        """
        scripts/layer_go.py --classification_file {input.classification_file} --output {output} --biomart_host {params.biomart_host} --biomart_dataset {params.biomart_dataset} --db {params.db} &> {log}
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
        os.path.join(path_output, "logs", prefix, "{db}", "layer_interproscan.log")
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
        os.path.join(path_output, "logs", prefix, "{db}", "layer_exons.log")
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
        os.path.join(path_output, "logs", prefix, "{db}", "layer_junctions.log")
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
        os.path.join(path_output, "logs", prefix, "{db}", "layer_nmd.log")
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
        species=config["species"],
        db=config.get("db")
    output:
        _output_layer_db("layer_reactome", external_rule=rules.get_reactome.output)
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "layer_reactome.log")
    shell:
        """
        scripts/layer_reactome.py --reactome_file {input.reactome_file} --classification_file {input.classification_file} --biomart_host {params.biomart_host} --biomart_dataset {params.biomart_dataset} --species {params.species:q} --db {params.db} --output {output} &> {log}
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
        os.path.join(path_output, "logs", prefix, "{db}", "layer_repeatmasker.log")
    shell:
        """
        scripts/layer_repeatmasker.py --keep_version {params.keep_version} --repeatmasker_file {input.repeatmasker_file} --classification_file {input.classification_file} --output {output} &> {log}
        """

rule layer_mirna_bs:
    input:
        mirna_bs_file=rules.get_mirna_bs_annotation.output
    output:
        _output_layer_db("layer_mirna_bs", external_rule=rules.get_mirna_bs_annotation.output)
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "layer_mirna_bs.log")
    shell:
        """
        sort -u -V -k1,1 -k4,4n -k5,5n {input.mirna_bs_file} > {output} 2> {log}
        """

rule layer_uniprot:
    input:
        uniprot_file=rules.get_uniprot_phosphosite_annotation.output,
    output:
        _output_layer_db("layer_uniprot", external_rule=rules.get_uniprot_phosphosite_annotation.output)
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "layer_uniprot.log")
    shell:
        """
        sort -u -V -k1,1 -k3,3 -k4,4n -k5,5n {input.uniprot_file} > {output} 2> {log}
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
        os.path.join(path_output, "logs", prefix, "{db}", "layer_utrscan.log")
    shell:
        """
        scripts/layer_utrscan.py --keep_version {params.keep_version} --utrscan_file {input.utrscan_file} --classification_file {input.classification_file} --output {output} &> {log}
        """

rule layer_nls:
    input:
        nls = rules.expand_nls.output,
        classification=select_sqanti_classification
    output: 
       _output_layer_db("layer_nls", external_rule=rules.parse_nls.output)
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "layer_nls.log")
    shell:
        """
        scripts/layer_nls.py \
            --nls_file {input.nls} \
            --classification_file {input.classification} \
            --output {output}.tmp &> {log}
        sort -V -k1,1 -k4,4n {output}.tmp > {output}
        rm {output}.tmp
        """

# GET EVERYTHING TOGETHER
rule tappas_annotation:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        transcript_block = [
            rules.layer_utrscan.output,
            rules.layer_repeatmasker.output if is_feature_enabled("repeat_masker", is_feature_enabled("repeatmasker", True)) else [],
            rules.layer_nmd.output,
            rules.layer_mirna_bs.output if config.get("mirna_regions_to_use") else []
        ] + config.get("transcript_gtf", []),
        genomic_block = [
            rules.layer_exons.output,
            rules.layer_junctions.output,
        ] + config.get("genomic_gtf", []),
        protein_block = [
            rules.layer_go.output if config.get("layer_go", "no") == "si" else [],
            rules.layer_reactome.output if config.get("reactome") else [],
            rules.layer_interproscan.output,
            rules.layer_uniprot.output,
            rules.layer_nls.output if config.get("nls_model") else []
        ] + config.get("protein_gtf", []),
        classification_file=select_sqanti_classification,
        gene_desc=[],
        protein_assoc=select_prot_assoc
    output:
        os.path.join(path_output, "data", prefix, f"{species_name}_tappas_{{db}}_annotation_file.gff3")
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "tappas_annotation.log")
    shell:
        """
        scripts/t2goAnnotationFile.py --classification_file {input.classification_file}  \
         --gene_desc_file {input.gene_desc} --input_transcripts {input.transcript_block} --input_genomic {input.genomic_block} \
         --input_protein  {input.protein_block} --output {output}.tmp --gene_desc_file {input.gene_desc} --protein_association {input.protein_assoc} &> {log}
        sort -V -k1,1 -k4,4n {output}.tmp > {output}
        rm {output}.tmp
        """


rule renameFeatures:
    conda:
        "../../envs/isoannotpy.yaml"
    input:
        rules.tappas_annotation.output
    output:
        os.path.join(path_output, "data", prefix, f"{species_name}_tappas_{{db}}_annotation_file.gff3_mod") 
    log:
        os.path.join(path_output, "logs", prefix, "{db}", "renameFeatures.log")
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
