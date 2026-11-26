import gzip
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

from transcript2reference import group_exons_by_transcript, _retrieve_pacbio_info


def exon_line(transcript_id, start, end, strand="+", chrom="1"):
    return f'{chrom}\tPacBio\texon\t{start}\t{end}\t.\t{strand}\t.\ttranscript_id "{transcript_id}"; gene_id "PB.1";\n'


# PB.1.1 has its exons apart in the file, with those of PB.2.1 in between
GTF_LINES = [exon_line("PB.1.1", 101, 200), exon_line("PB.2.1", 1101, 1200, "-"), exon_line("PB.1.1", 301, 400),
             exon_line("PB.2.1", 1301, 1400, "-"), exon_line("PB.1.1", 501, 600)]


def test_group_exons_by_transcript():
    """Every exon of a transcript is found, also when other transcripts' exons are in between, in file order."""
    assert group_exons_by_transcript(GTF_LINES) == {
        "PB.1.1": [GTF_LINES[0], GTF_LINES[2], GTF_LINES[4]],
        "PB.2.1": [GTF_LINES[1], GTF_LINES[3]],
    }


def test_retrieve_pacbio_info():
    """
    CDS and NMD of each isoform, the same as the per-isoform scan of the whole GTF gave.
    PB.1.1: CDS 111-200 in its first exon, stop more than 55 nt before the last junction, so NMD.
    PB.2.1: minus strand, CDS from 1380 down to 1200, across the junction into the last exon, so not NMD.
    """
    orfs = {"PB.1.1\tPB.1.1.p1|30_aa|0.9|complete|11|100": "M" * 29 + "*",
            "PB.2.1\tPB.2.1.p1|27_aa|0.9|complete|21|101": "M" * 26 + "*"}
    assert _retrieve_pacbio_info(GTF_LINES, orfs) == {
        "PB.1.1": {"start": [111], "end": [200], "len": 90, "nmd": "NMD", "chr": "1"},
        "PB.2.1": {"start": [1200, 1301], "end": [1200, 1380], "len": 81, "nmd": "nonNMD", "chr": "1"},
    }


@pytest.mark.parametrize("name", ["isoforms.fna", "isoforms.fas", "isoforms.fna.gz"])
def test_check_chromosomes_fasta_any_name(run_script, tmp_path, name):
    """A FASTA is recognised by its content: .fna or .fas were read as a GTF and stopped the run with a naming error."""
    user_input = tmp_path / name
    opener = gzip.open if name.endswith(".gz") else open
    with opener(user_input, "wt") as h:
        h.write(">PB.1.1\nACGTACGT\n")
    result = run_script("check_chromosomes.py", "--user_input", user_input,
                        "--ref_gtf", tmp_path / "missing.gtf", "--ref_fasta", tmp_path / "missing.fa")
    assert result.returncode == 0, result.stderr
    assert result.stderr == ""


def test_check_chromosomes_gtf_naming_error(run_script, tmp_path):
    """A GTF with chr1 against a genome with 1 still stops with the naming error."""
    (tmp_path / "user.gtf").write_text(exon_line("PB.1.1", 101, 200, chrom="chr1"))
    (tmp_path / "ref.gtf").write_text(exon_line("ENST1", 101, 200))
    (tmp_path / "genome.fa").write_text(">1 dna:chromosome\nACGT\n")
    result = run_script("check_chromosomes.py", "--user_input", tmp_path / "user.gtf",
                        "--ref_gtf", tmp_path / "ref.gtf", "--ref_fasta", tmp_path / "genome.fa")
    assert result.returncode == 1
    assert "Transcriptome IDs: ['chr1'] ... NOT found in Reference FASTA: ['1']" in result.stderr


def test_check_chromosome_consistency_log(tmp_path):
    """The check writes its messages to a log, like the other rules, instead of only to the Snakemake output."""
    if not SNAKEMAKE or not os.path.isfile(SNAKEMAKE):
        pytest.skip("Snakemake not found (set ISOANNOT_SNAKEMAKE or put snakemake on PATH)")
    ref_dir = tmp_path / "data" / "Hsapiens" / "config" / "ensembl"
    ref_dir.mkdir(parents=True)
    (ref_dir / "ref.gtf").write_text(exon_line("ENST1", 101, 200))
    (ref_dir / "genome.fa").write_text(">1\nACGT\n>2\nACGT\n")
    user_gtf = tmp_path / "isoforms.gtf"
    user_gtf.write_text(exon_line("PB.1.1", 101, 200) + exon_line("PB.2.1", 101, 200, chrom="2"))
    config_dir = os.path.join(ISOANNOT_DIR, "config", "mytranscripts", "hsapiens")
    cmd = [SNAKEMAKE, "-p", "--snakefile", os.path.join(config_dir, "Snakefile.smk"),
           "--configfile", os.path.join(config_dir, "config.yaml"),
           "--config", "db=mytranscripts", f"path_output={tmp_path}", f"gtf_cdna={user_gtf}",
           "ensembl_gtf=https://host/p/ref.gtf.gz", "ensembl_reference=https://host/p/genome.fa.gz",
           "--directory", ISOANNOT_DIR, "--cores", "1", "--allowed-rules", "check_chromosome_consistency",
           "--nolock", "--drop-metadata", str(tmp_path / "data" / "Hsapiens" / "config" / "check_chromosomes.done")]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    log = (tmp_path / "logs" / "Hsapiens" / "check_chromosomes.log").read_text()
    assert "WARNING: Only 50.00% of your transcriptome IDs are in the reference GTF." in log
