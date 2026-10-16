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
