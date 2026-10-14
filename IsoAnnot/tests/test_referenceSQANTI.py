import csv
import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_reference_sqanti(run_script, tmp_path, gtf, genome, chr2acc, database):
    outputs = {name: tmp_path / f"{database}_sqanti_{name}.txt" for name in ("classification", "junctions", "NMD")}
    result = run_script("referenceSQANTI.py", "--gtf_file", gtf, "--reference_file", genome, "--chr_ref", chr2acc,
                        "--database", database, "--output_classification", outputs["classification"],
                        "--output_junctions", outputs["junctions"], "--output_nmd", outputs["NMD"])
    return result, outputs


def test_refseq_mode_writes_every_transcript(run_script, synthetic_dir, tmp_path):
    """RefSeq GTF on NC_000001.11 with the Ensembl genome "1": all six transcripts are written on "1"."""
    result, outputs = run_reference_sqanti(run_script, tmp_path, synthetic_dir / "refseq.gtf",
                                           synthetic_dir / "genome.fa", synthetic_dir / "chr2acc", "refseq")
    assert result.returncode == 0, result.stderr
    rows = read_tsv(outputs["classification"])
    assert [(r["isoform"], r["chrom"], r["strand"], r["length"], r["exons"], r["associated_gene"],
             r["CDS_start"], r["CDS_end"], r["refProt"]) for r in rows] == [
        ("NM_900001.1", "1", "+", "750", "3", "SYNA", "101", "550", "NP_900001.1"),
        ("NM_900002.1", "1", "+", "850", "4", "SYNA", "101", "361", "NP_900002.1"),
        ("NM_900003.2", "1", "-", "710", "4", "SYNB", "51", "611", "NP_900003.1"),
        ("NM_900004.1", "1", "-", "770", "5", "SYNB", "111", "671", "NP_900004.1"),
        ("NR_900005.1", "1", "+", "700", "2", "SYNC", "NA", "NA", "NA"),
        ("NM_900006.3", "1", "+", "1100", "1", "SYND", "101", "901", "NP_900006.2"),
    ]
    junctions = read_tsv(outputs["junctions"])
    assert len(junctions) == 13
    assert {(j["chrom"], j["splice_site"], j["canonical"]) for j in junctions} == {("1", "GTAG", "canonical")}
    assert outputs["NMD"].read_text().splitlines() == [
        "NM_900001.1\tnonNMD", "NM_900002.1\tNMD", "NM_900003.2\tnonNMD",
        "NM_900004.1\tnonNMD", "NR_900005.1\tNA", "NM_900006.3\tnonNMD"]


def test_refseq_mode_matches_the_ensembl_mode(run_script, synthetic_dir, tmp_path):
    """Same models, so the refseq and ensembl classifications only differ in IDs and biotype."""
    _, outputs = run_reference_sqanti(run_script, tmp_path, synthetic_dir / "refseq.gtf",
                                      synthetic_dir / "genome.fa", synthetic_dir / "chr2acc", "refseq")
    ensembl = read_tsv(synthetic_dir / "ensembl_sqanti_classification.txt")
    columns = ["chrom", "strand", "length", "exons", "coding", "ORF_length", "CDS_length", "CDS_start", "CDS_end"]
    assert [[r[c] for c in columns] for r in read_tsv(outputs["classification"])] == \
        [[r[c] for c in columns] for r in ensembl]


@pytest.mark.parametrize("database, lookups", [("ensembl", 66), ("refseq", 54)])
def test_genome_with_other_chromosome_names_stops_the_run(run_script, synthetic_dir, tmp_path, database, lookups):
    """A genome on "chr1" matches no transcript: the run fails instead of writing empty tables with exit 0."""
    genome = tmp_path / "genome_chr.fa"
    genome.write_text((synthetic_dir / "genome.fa").read_text().replace(">1 ", ">chr1 ", 1))
    result, outputs = run_reference_sqanti(run_script, tmp_path, synthetic_dir / f"{database}.gtf",
                                           genome, synthetic_dir / "chr2acc", database)
    assert result.returncode != 0
    # Every GTF line the parser looks up fails: 66 non-gene lines in ensembl.gtf, all 54 lines in refseq.gtf
    assert f"{lookups}/{lookups} chromosome lookups failed" in result.stderr
    assert [path for path in outputs.values() if path.exists()] == []
