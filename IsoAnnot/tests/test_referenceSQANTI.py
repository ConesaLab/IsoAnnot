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


def test_ambiguous_splice_site_bases_do_not_truncate_the_junctions(run_script, tmp_path):
    """An R or Y at a splice site is read as N and processing goes on; it used to stop with exit 0."""
    genome = ["A"] * 120
    for pos, base in {11: "G", 12: "R", 19: "A", 20: "G",    # TA, + strand: GR..AG
                      51: "G", 52: "T", 59: "A", 60: "G",    # TB, + strand: GT..AG
                      91: "C", 92: "Y", 99: "A", 100: "C"}.items():  # TC, - strand: CY..AC
        genome[pos - 1] = base
    genome_file = tmp_path / "genome.fa"
    genome_file.write_text(">1\n" + "".join(genome) + "\n")
    gtf = tmp_path / "ensembl.gtf"
    lines = []
    for transcript, strand, exons in [("TA", "+", [(1, 10), (21, 30)]), ("TB", "+", [(41, 50), (61, 70)]),
                                      ("TC", "-", [(81, 90), (101, 110)])]:
        attributes = f'gene_id "G{transcript}"; transcript_id "{transcript}"; gene_name "G{transcript}"; transcript_biotype "lncRNA";'
        lines.append(f"1\ttest\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\t{attributes}")
        lines.extend(f"1\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\t{attributes}" for start, end in exons)
    gtf.write_text("\n".join(lines) + "\n")
    chr2acc = tmp_path / "chr2acc"
    chr2acc.write_text("#Chromosome\tAccession.version\n1\tNC_000001.11\n")

    result, outputs = run_reference_sqanti(run_script, tmp_path, gtf, genome_file, chr2acc, "ensembl")
    assert result.returncode == 0, result.stderr
    assert [(j["isoform"], j["splice_site"], j["canonical"]) for j in read_tsv(outputs["junctions"])] == [
        ("TA", "GNAG", "non_canonical"), ("TB", "GTAG", "canonical"), ("TC", "GTNG", "non_canonical")]
    assert result.stderr.count("2 splice sites have bases other than A, C, G, T or N") == 1


def test_stopped_run_leaves_no_truncated_junctions_file(synthetic_dir, tmp_path, monkeypatch):
    """A run that stops in the middle of the junctions leaves no junctions file, not a header-only one."""
    import referenceSQANTI

    def stop(site):
        raise RuntimeError("stopped at a - strand junction")
    # The first - strand junction comes after the + strand junctions of SYNA have been written
    monkeypatch.setattr(referenceSQANTI, "ReverseComplement", stop)
    junctions = tmp_path / "ensembl_sqanti_junctions.txt"
    monkeypatch.setattr(sys, "argv", ["referenceSQANTI.py", "--gtf_file", str(synthetic_dir / "ensembl.gtf"),
                                      "--reference_file", str(synthetic_dir / "genome.fa"),
                                      "--chr_ref", str(synthetic_dir / "chr2acc"), "--database", "ensembl",
                                      "--output_classification", str(tmp_path / "classification.txt"),
                                      "--output_junctions", str(junctions), "--output_nmd", str(tmp_path / "NMD.txt")])
    with pytest.raises(RuntimeError, match="stopped at a - strand junction"):
        referenceSQANTI.main()
    assert not junctions.exists()


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
