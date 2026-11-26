import os
import sys

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)


def test_refseq_exons_get_the_ensembl_chromosome(run_script, synthetic_dir, tmp_path):
    """RefSeq exons on NC_000001.11 get Chr=chr1 (regression guard for the target_db rename)."""
    output = tmp_path / "layer_exons.gtf"
    result = run_script("layer_exons.py", "--gtf_file", synthetic_dir / "refseq.gtf",
                        "--chr_ref", synthetic_dir / "chr2acc", "--output", output)
    assert result.returncode == 0, result.stderr
    lines = output.read_text().splitlines()
    assert len(lines) == 19
    assert {line.split("\t")[8] for line in lines} == {"Chr=chr1"}
    assert lines[0] == "NM_900001.1\ttappAS\texon\t1001\t1200\t.\t+\t.\tChr=chr1"


def test_unmapped_chromosomes_stop_the_layer(run_script, synthetic_dir, tmp_path):
    """The failure check runs after the lookups, so a chr2acc without NC_000001.11 fails with no output."""
    chr2acc = tmp_path / "chr2acc"
    chr2acc.write_text("#Chromosome\tAccession.version\n2\tNC_000002.12\n")
    output = tmp_path / "layer_exons.gtf"
    result = run_script("layer_exons.py", "--gtf_file", synthetic_dir / "refseq.gtf",
                        "--chr_ref", chr2acc, "--output", output)
    assert result.returncode != 0
    assert "19/19 chromosome lookups failed" in result.stderr
    assert not output.exists()
    assert not (tmp_path / "layer_exons.gtf.tmp").exists()
