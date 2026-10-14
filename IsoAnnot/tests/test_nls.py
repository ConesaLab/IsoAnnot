import os
import subprocess
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

ISOANNOT_DIR = os.path.dirname(main_path)
SNAKEMAKE = os.environ.get("ISOANNOT_SNAKEMAKE", "/media/pabloati/HDD1/conda/envs/snakemake/bin/snakemake")

# First two lines of every NucImport output, as NucImportMay2012.jar writes them
NUCIMPORT_HEADER = ("Protein ID\tImport Probability\t cNLS probability\tPredicted Most probable cNLS\t"
                    "Predicted most probable cNLS class number\tPosition")
NUCIMPORT_STARS = "*" * 93


@pytest.fixture
def nls_dir(tmp_path):
    """nls folder of the ensembl hsapiens run whose path_output is tmp_path."""
    return tmp_path / "data" / "Hsapiens" / "output" / "ensembl" / "nls"


@pytest.fixture
def run_rule(tmp_path):
    """
    Runs the Snakefile rules for the given targets, with the ensembl hsapiens config and path_output=tmp_path.
    Only the given rule is allowed, so its input files are used as they are and nothing upstream runs.
    Scripts run with the test interpreter (first python3 on PATH). Returns the CompletedProcess.
    """
    if not os.path.isfile(SNAKEMAKE):
        pytest.skip(f"Snakemake not found at {SNAKEMAKE} (set ISOANNOT_SNAKEMAKE to use another one)")

    def _run(rule, targets, *config, dry_run=False):
        config_dir = os.path.join(ISOANNOT_DIR, "config", "ensembl", "hsapiens")
        cmd = [SNAKEMAKE, "-p",
               "--snakefile", os.path.join(config_dir, "Snakefile.smk"),
               "--configfile", os.path.join(config_dir, "config.yaml"),
               "--config", "db=ensembl", f"path_output={tmp_path}", *config,
               "--directory", ISOANNOT_DIR, "--cores", "2", "--allowed-rules", rule,
               "--nolock", "--drop-metadata", *[str(t) for t in targets]]
        if dry_run:
            cmd.insert(1, "-n")
        env = dict(os.environ, PATH=os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"])
        return subprocess.run(cmd, capture_output=True, text=True, env=env)
    return _run


def nucimport_row(prot_id):
    return f"{prot_id}\t0.9949979710297792\t0.806466703423478\tAALLLLKRKRRSPPKKKRKV\tClass1\t15\t"


def nucimport_output(*prot_ids):
    """NucImport output: a header line, a line of stars, then one line per protein."""
    return "\n".join([NUCIMPORT_HEADER, NUCIMPORT_STARS] + [nucimport_row(p) for p in prot_ids]) + "\n"


def test_merge_nls_chunks_keeps_first_protein(run_rule, nls_dir):
    """The merge skipped 3 lines of every chunk after the first, so it lost their first protein."""
    chunk_dir = nls_dir / "tmp"
    chunk_dir.mkdir(parents=True)
    (chunk_dir / "output_chunk_0.txt").write_text(nucimport_output("P1", "P4"))
    (chunk_dir / "output_chunk_1.txt").write_text(nucimport_output("P2", "P5"))
    (chunk_dir / "output_chunk_2.txt").write_text(nucimport_output("P3"))
    merged = nls_dir / "output_merged.txt"
    result = run_rule("merge_nls_chunks", [merged], "nls_chunks=3")
    assert result.returncode == 0, result.stdout + result.stderr
    assert merged.read_text() == nucimport_output("P1", "P4", "P2", "P5", "P3")
