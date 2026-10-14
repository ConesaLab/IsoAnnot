import os
import subprocess
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

ISOANNOT_DIR = os.path.dirname(main_path)
SNAKEMAKE = os.environ.get("ISOANNOT_SNAKEMAKE", "/media/pabloati/HDD1/conda/envs/snakemake/bin/snakemake")


@pytest.fixture
def dry_run(tmp_path):
    """
    Runs snakemake -n with the arguments isoannot.sh builds, for the hsapiens config of a database.
    --use-conda is left out because a dry run does not need the environments. The output folder is a
    fresh tmp_path. Returns the CompletedProcess, with the planned commands (-p) in stdout.
    """
    if not os.path.isfile(SNAKEMAKE):
        pytest.skip(f"Snakemake not found at {SNAKEMAKE} (set ISOANNOT_SNAKEMAKE to use another one)")

    def _run(database, *config, configfile=None):
        config_dir = os.path.join(ISOANNOT_DIR, "config", database, "hsapiens")
        cmd = [SNAKEMAKE, "-n", "-p",
               "--snakefile", os.path.join(config_dir, "Snakefile.smk"),
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
    assert f"zcat -f {refseq_dir}/GCF_1_genomic.gtf|grep -v -F 'partial=true' | grep '^NC_' > {refseq_dir}/GCF_1_genomic_nopartial_nc.gtf" in planned


def test_dry_run_mytranscripts_fasta(dry_run, synthetic_dir):
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


@pytest.mark.parametrize("key, layer", [
    ("layer_go", "layer_go"),
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
