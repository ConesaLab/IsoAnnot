import os
import shutil
import subprocess
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

ISOANNOT_DIR = os.path.dirname(main_path)
# ISOANNOT_SNAKEMAKE, or snakemake on PATH
SNAKEMAKE = os.environ.get("ISOANNOT_SNAKEMAKE") or shutil.which("snakemake")


@pytest.fixture
def dry_run(tmp_path):
    """
    Runs snakemake -n with the arguments isoannot.sh builds, for the hsapiens config of a database.
    --use-conda is left out because a dry run does not need the environments. The output folder is a
    fresh tmp_path. Returns the CompletedProcess, with the planned commands (-p) in stdout.
    """
    if not SNAKEMAKE or not os.path.isfile(SNAKEMAKE):
        pytest.skip("Snakemake not found (set ISOANNOT_SNAKEMAKE or put snakemake on PATH)")

    def _run(database, *config, snakefile=None, configfile=None):
        config_dir = os.path.join(ISOANNOT_DIR, "config", database, "hsapiens")
        cmd = [SNAKEMAKE, "-n", "-p",
               "--snakefile", snakefile or os.path.join(config_dir, "Snakefile.smk"),
               "--configfile", configfile or os.path.join(config_dir, "config.yaml"),
               "--config", f"db={database}", f"path_output={tmp_path}", *config,
               "--directory", ISOANNOT_DIR, "--cores", "8", "all",
               "--rerun-incomplete", "--nolock", "--resources", "n_downloads=1"]
        return subprocess.run(cmd, capture_output=True, text=True)
    return _run


def job_counts(output):
    """Rule name -> number of jobs, from the first "Job stats" table of a dry run."""
    lines = output.splitlines()
    start = lines.index("Job stats:") + 3
    counts = {}
    for line in lines[start:]:
        if not line.strip():
            break
        rule, count = line.split()
        counts[rule] = int(count)
    return counts


@pytest.mark.parametrize("database", ["ensembl", "refseq"])
def test_dry_run_reference_modes(dry_run, database):
    """Reference modes build the classification with referenceSQANTI and reach the final GFF3."""
    result = dry_run(database)
    assert result.returncode == 0, result.stdout + result.stderr
    jobs = job_counts(result.stdout)
    assert jobs["run_refsqanti"] == 1
    assert jobs["tappas_annotation"] == 1
    assert jobs["renameFeatures"] == 1
    assert "run_sqanti" not in jobs
    assert "transcript_to_reference" not in jobs
    assert jobs[f"prepare_{database}_gtf"] == 1


@pytest.mark.parametrize("database", ["ensembl", "refseq"])
def test_dry_run_plain_urls(dry_run, tmp_path, database):
    """URLs without .gz are downloaded and used as they are, instead of crashing at parse time."""
    result = dry_run(database,
                     "ensembl_proteins=https://host/p/Homo_sapiens.GRCh38.pep.all.fa",
                     "ensembl_cdna=https://host/p/Homo_sapiens.GRCh38.cdna.all.fa",
                     "ensembl_gtf=https://host/p/Homo_sapiens.GRCh38.115.chr.gtf",
                     "refseq_proteins=https://host/p/GCF_1_protein.faa",
                     "refseq_cdna=https://host/p/GCF_1_rna.fna",
                     "refseq_gtf=https://host/p/GCF_1_genomic.gtf")
    assert result.returncode == 0, result.stdout + result.stderr
    jobs = job_counts(result.stdout)
    assert "prepare_ensembl_gtf" not in jobs
    assert jobs["get_ensembl_gtf"] == 1
    assert jobs["get_ensembl_proteins"] == 1
    assert jobs["prepare_refseq_gtf"] == 1
    planned = " ".join(result.stdout.split())
    ensembl_dir = tmp_path / "data" / "Hsapiens" / "config" / "ensembl"
    refseq_dir = tmp_path / "data" / "Hsapiens" / "config" / "refseq"
    assert f"gunzip -cf {ensembl_dir}/Homo_sapiens.GRCh38.pep.all.fa > {ensembl_dir}/ensembl_proteins.fa.tmp" in planned
    assert f"zcat -f {refseq_dir}/GCF_1_genomic.gtf | grep -v -F 'partial=true' | awk " in planned


def test_dry_run_plain_uniprot_urls(dry_run, tmp_path, synthetic_dir):
    """UniProt URLs without .gz are passed to the scripts, instead of being dropped by a .dat.gz/.fasta.gz filter."""
    configfile = tmp_path / "plain_uniprot.yaml"
    configfile.write_text("uniprot_fasta:\n  - https://host/p/UP1_9606.fasta\n"
                          "uniprot_dat:\n  - https://host/p/UP1_9606.dat\n")
    # mytranscripts, because parse_protein_databases only runs in that mode
    result = dry_run("mytranscripts", f"fasta_cdna={synthetic_dir / 'ensembl_cdna.fa'}", configfile=str(configfile))
    assert result.returncode == 0, result.stdout + result.stderr
    uniprot_dir = tmp_path / "data" / "Hsapiens" / "config" / "uniprot"
    assert f"--uniprot_files {uniprot_dir}/UP1_9606.dat --output " in planned_command(result.stdout, "uniprot_parse.py")
    assert f"--uniprot_fasta {uniprot_dir}/UP1_9606.fasta --ensembl_fasta " in planned_command(result.stdout, "parse_protein_databases.py")
    assert f"--uniprot_fasta {uniprot_dir}/UP1_9606.fasta --refseq_fasta " in planned_command(result.stdout, "uniprotPhosphosite_genomicCoordinates.py")


def planned_command(output, script):
    """The planned shell command (-p) that runs scripts/<script>, as one line."""
    blocks = output.split("\n\n")
    return next(" ".join(block.split()) for block in blocks if f"scripts/{script}" in block)


@pytest.mark.parametrize("database", ["ensembl", "refseq"])
def test_dry_run_reference_protein_ids(dry_run, database, tmp_path):
    """Reference modes give t2goAnnotationFile.py the same-sequence protein IDs (GitHub #10)."""
    result = dry_run(database)
    assert result.returncode == 0, result.stdout + result.stderr
    jobs = job_counts(result.stdout)
    assert jobs["reference_protein_assoc"] == 1
    assert jobs["parse_protein_databases"] == 1
    assoc = tmp_path / "data" / "Hsapiens" / "output" / database / "reference_protein_assoc_data.txt"
    assert f"--protein_association {assoc}" in planned_command(result.stdout, "t2goAnnotationFile.py")
    assert f"--database {database} --output {assoc}" in planned_command(result.stdout, "reference_protein_assoc.py")


def test_dry_run_tappas_annotation_keeps_line_order(dry_run, tmp_path):
    """
    The GFF3 keeps the order written by t2goAnnotationFile.py. sort -V -k1,1 -k4,4n moved GO, gene and
    protein lines before the transcript line. tappAS needs the transcript line first in each block.
    """
    result = dry_run("ensembl")
    assert result.returncode == 0, result.stdout + result.stderr
    command = planned_command(result.stdout, "t2goAnnotationFile.py")
    gff3 = tmp_path / "data" / "Hsapiens" / "human_tappas_ensembl_annotation_file.gff3"
    assert command.endswith(f"mv {gff3}.tmp {gff3}")
    assert "sort" not in command


def test_dry_run_mytranscripts_fasta(dry_run, synthetic_dir, tmp_path):
    """mytranscripts with a FASTA input runs SQANTI3 with a GMAP index and transcript2reference."""
    cdna = synthetic_dir / "ensembl_cdna.fa"
    result = dry_run("mytranscripts", f"fasta_cdna={cdna}")
    assert result.returncode == 0, result.stdout + result.stderr
    jobs = job_counts(result.stdout)
    assert jobs["run_sqanti"] == 1
    assert jobs["run_gmap_index"] == 1
    assert jobs["transcript_to_reference"] == 1
    assert jobs["tappas_annotation"] == 1
    assert "run_refsqanti" not in jobs
    assert f"--isoforms {cdna}" in result.stdout
    # Protein IDs come from transcript2reference.py, not reference_protein_assoc
    assert "reference_protein_assoc" not in jobs
    assoc = tmp_path / "data" / "Hsapiens" / "output" / "mytranscripts" / "protein_assoc_data.txt"
    assert f"--protein_association {assoc}" in planned_command(result.stdout, "t2goAnnotationFile.py")


@pytest.mark.parametrize("key, layer", [
    ("layer_go", "layer_go"),
    ("reactome", "layer_reactome"),
])
def test_dry_run_empty_toggle_is_off(dry_run, tmp_path, key, layer):
    """An empty value in the species YAML keeps the layer off, as it was before is_feature_enabled."""
    configfile = tmp_path / "empty_toggle.yaml"
    configfile.write_text(f'{key}: ""\n')
    result = dry_run("ensembl", configfile=str(configfile))
    assert result.returncode == 0, result.stdout + result.stderr
    jobs = job_counts(result.stdout)
    assert layer not in jobs
    assert jobs["tappas_annotation"] == 1


@pytest.mark.parametrize("value, enabled", [
    ("si", True),
    ("True", True),
    ("true", True),
    ("yes", True),
    ("no", False),
    ("False", False),
])
def test_dry_run_layer_go_toggle(dry_run, value, enabled):
    """layer_go accepts si/true/yes to turn the layer on and no/false to turn it off."""
    result = dry_run("ensembl", f"layer_go={value}")
    assert result.returncode == 0, result.stdout + result.stderr
    jobs = job_counts(result.stdout)
    assert ("layer_go" in jobs) == enabled
    assert jobs["tappas_annotation"] == 1


@pytest.mark.parametrize("value, enabled", [
    ("https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt", True),
    ("no", False),
    ("False", False),
])
def test_dry_run_reactome_toggle(dry_run, value, enabled):
    """reactome=no (or False, what YAML makes of "reactome: no") turns the Reactome layer off, as the BioMart error advises."""
    result = dry_run("ensembl", f"reactome={value}")
    assert result.returncode == 0, result.stdout + result.stderr
    jobs = job_counts(result.stdout)
    assert ("layer_reactome" in jobs) == enabled
    assert ("get_reactome" in jobs) == enabled


@pytest.mark.parametrize("config, enabled", [
    ([], True),
    (["repeat_masker=False"], False),
])
def test_dry_run_repeat_masker_toggle(dry_run, config, enabled):
    """RepeatMasker is on unless repeat_masker=False, and Dfam is only downloaded when it is on."""
    result = dry_run("ensembl", *config)
    assert result.returncode == 0, result.stdout + result.stderr
    jobs = job_counts(result.stdout)
    assert ("run_repeatmasker" in jobs) == enabled
    assert ("layer_repeatmasker" in jobs) == enabled
    assert ("get_dfam" in jobs) == enabled
    assert ("install_repeatmasker" in jobs) == enabled
    assert jobs["tappas_annotation"] == 1


@pytest.mark.parametrize("config, dfam_files, uncurated", [
    ([], ["dfam40.0.h5", "dfam40.curated.consensus.0.h5"], False),
    (["repeatmasker_dfam=all"], ["dfam40.0.h5", "dfam40.curated.consensus.0.h5",
                                 "dfam40.uncurated.consensus.0.h5", "dfam40.uncurated.consensus.1.h5"], True),
])
def test_dry_run_repeatmasker_dfam(dry_run, tmp_path, config, dfam_files, uncurated):
    """The curated Dfam files by default; repeatmasker_dfam=all adds the uncurated ones and RepeatMasker -uncurated."""
    result = dry_run("ensembl", "repeat_masker=True", *config)
    assert result.returncode == 0, result.stdout + result.stderr
    jobs = job_counts(result.stdout)
    assert jobs["get_dfam"] == len(dfam_files)
    for dfam_file in dfam_files:
        assert f"{tmp_path}/data/global/dfam/{dfam_file}" in result.stdout
    # The command is split over two lines in the rule
    lines = result.stdout.splitlines()
    n = next(i for i, line in enumerate(lines) if line.strip().startswith(f"{tmp_path}/data/global/repeatmasker/RepeatMasker "))
    repeatmasker_command = lines[n] + lines[n + 1]
    assert ("-uncurated" in repeatmasker_command) == uncurated
    assert f"-famdb_data_dir {tmp_path}/data/global/dfam " in repeatmasker_command
    assert f"-cache_dir {tmp_path}/data/global/dfam/cache " in repeatmasker_command
    assert "-pa 2" in repeatmasker_command  # 8 cores in the dry run, 4 per RMBlast job


def test_dry_run_repeatmasker_dfam_invalid(dry_run):
    """A repeatmasker_dfam value other than curated or all stops the run and names the valid values."""
    result = dry_run("ensembl", "repeat_masker=True", "repeatmasker_dfam=uncurated")
    assert result.returncode != 0
    assert "repeatmasker_dfam must be curated or all, not 'uncurated'." in result.stdout + result.stderr


@pytest.mark.parametrize("config, score_option", [
    ([], "--score 0.95 "),
    (["mirna_db_evidence_score_threshold="], ""),
])
def test_dry_run_mirna_score_threshold(dry_run, config, score_option):
    """An empty mirna_db_evidence_score_threshold leaves the --score option out instead of giving it no value."""
    result = dry_run("ensembl", *config)
    assert result.returncode == 0, result.stdout + result.stderr
    command = planned_command(result.stdout, "filter_mirna_bs.py")
    assert f"--species 'Homo sapiens' {score_option}--mirwalk_output " in command


def test_dry_run_mirna_score_threshold_missing(dry_run, tmp_path):
    """A species config without mirna_db_evidence_score_threshold gave --score None."""
    species_config = os.path.join(ISOANNOT_DIR, "config", "ensembl", "hsapiens", "config.yaml")
    with open(species_config) as h:
        lines = [line for line in h if not line.startswith("mirna_db_evidence_score_threshold:")]
    configfile = tmp_path / "config.yaml"
    configfile.write_text("".join(lines))
    # The species Snakefile loads its own config.yaml, which has the key, so include the generic one directly
    snakefile = tmp_path / "Snakefile.smk"
    snakefile.write_text(f'configfile: "{configfile}"\n\n'
                         f'include: "{os.path.join(ISOANNOT_DIR, "config", "generic", "Snakefile_ensembl.smk")}"\n')
    result = dry_run("ensembl", snakefile=str(snakefile), configfile=str(configfile))
    assert result.returncode == 0, result.stdout + result.stderr
    command = planned_command(result.stdout, "filter_mirna_bs.py")
    assert "--species 'Homo sapiens' --mirwalk_output " in command
    assert "None" not in command


def test_dry_run_reactome_missing(dry_run, tmp_path):
    """A species config without the reactome key has the Reactome layer off, instead of a KeyError."""
    species_config = os.path.join(ISOANNOT_DIR, "config", "ensembl", "hsapiens", "config.yaml")
    with open(species_config) as h:
        lines = [line for line in h if not line.startswith("reactome:")]
    configfile = tmp_path / "config.yaml"
    configfile.write_text("".join(lines))
    # The species Snakefile loads its own config.yaml, which has the key, so include the generic one directly
    snakefile = tmp_path / "Snakefile.smk"
    snakefile.write_text(f'configfile: "{configfile}"\n\n'
                         f'include: "{os.path.join(ISOANNOT_DIR, "config", "generic", "Snakefile_ensembl.smk")}"\n')
    result = dry_run("ensembl", snakefile=str(snakefile), configfile=str(configfile))
    assert result.returncode == 0, result.stdout + result.stderr
    jobs = job_counts(result.stdout)
    assert "layer_reactome" not in jobs
    assert "get_reactome" not in jobs
    assert jobs["tappas_annotation"] == 1


def test_dry_run_biomart_table(dry_run, tmp_path):
    """With biomart_table set, layer_go and layer_reactome get the saved table instead of querying BioMart."""
    table = tmp_path / "biomart_table.tsv"
    table.write_text("Gene stable ID\tGene name\tGO term accession\tGO term name\tGO domain\n")
    result = dry_run("ensembl", f"biomart_table={table}")
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"--db ensembl --biomart_table {table} " in planned_command(result.stdout, "layer_go.py")
    assert f"--db ensembl --biomart_table {table} --output " in planned_command(result.stdout, "layer_reactome.py")


def test_dry_run_without_biomart_table(dry_run):
    """Without biomart_table the layers query BioMart as before."""
    result = dry_run("ensembl")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "--biomart_table" not in planned_command(result.stdout, "layer_go.py")
    assert "--biomart_table" not in planned_command(result.stdout, "layer_reactome.py")
