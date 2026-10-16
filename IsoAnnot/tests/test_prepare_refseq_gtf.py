import gzip
import os
import subprocess

import pytest

ISOANNOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SNAKEMAKE = os.environ.get("ISOANNOT_SNAKEMAKE", "/media/pabloati/HDD1/conda/envs/snakemake/bin/snakemake")

# A fly-like chr2acc: chromosome arms are NT_ accessions, not NC_
CHR2ACC = "#Chromosome\tAccession.version\nX\tNC_004354.4\n2L\tNT_033779.5\n"
GTF_LINES = [
    "NC_004354.4\tRefSeq\ttranscript\t1\t100\t.\t+\t.\tgene_id \"geneX\"; transcript_id \"NM_1.1\";",
    "NT_033779.5\tRefSeq\ttranscript\t1\t100\t.\t+\t.\tgene_id \"gene2L\"; transcript_id \"NM_2.1\";",
    "NW_001845431.1\tRefSeq\ttranscript\t1\t100\t.\t+\t.\tgene_id \"geneScaffold\"; transcript_id \"NM_3.1\";",
    "NC_004354.4\tRefSeq\ttranscript\t1\t100\t.\t+\t.\tgene_id \"genePartial\"; transcript_id \"NM_4.1\"; partial \"true\"; partial=true",
]


def test_prepare_refseq_gtf_keeps_chromosomes_listed_in_chr2acc(tmp_path):
    """
    Chromosomes listed in chr2acc are kept even when they are not NC_ (the fruit-fly arms 2L/2R/3L/3R are NT_
    and were dropped, losing most fly genes). Unlisted scaffolds and partial genes are still dropped.
    """
    if not os.path.isfile(SNAKEMAKE):
        pytest.skip(f"Snakemake not found at {SNAKEMAKE} (set ISOANNOT_SNAKEMAKE to use another one)")
    config_dir = os.path.join(ISOANNOT_DIR, "config", "refseq", "hsapiens")
    refseq_dir = tmp_path / "data" / "Hsapiens" / "config" / "refseq"
    refseq_dir.mkdir(parents=True)
    with gzip.open(refseq_dir / "GCF_000001405.40_GRCh38.p14_genomic.gtf.gz", "wt") as handle:
        handle.write("#gtf-version 2.2\n" + "\n".join(GTF_LINES) + "\n")
    (refseq_dir / "chr2acc").write_text(CHR2ACC)
    target = refseq_dir / "GCF_000001405.40_GRCh38.p14_genomic_nopartial_nc.gtf"
    result = subprocess.run([SNAKEMAKE, "--snakefile", os.path.join(config_dir, "Snakefile.smk"),
                             "--configfile", os.path.join(config_dir, "config.yaml"),
                             "--config", "db=refseq", f"path_output={tmp_path}",
                             "--directory", ISOANNOT_DIR, "--cores", "1", "--nolock", str(target)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert target.read_text().splitlines() == GTF_LINES[:2]
