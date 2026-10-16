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

    def _run(database, *config):
        config_dir = os.path.join(ISOANNOT_DIR, "config", database, "hsapiens")
        cmd = [SNAKEMAKE, "-n", "-p",
               "--snakefile", os.path.join(config_dir, "Snakefile.smk"),
               "--configfile", os.path.join(config_dir, "config.yaml"),
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
