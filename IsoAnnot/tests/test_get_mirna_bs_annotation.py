import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

pytest.importorskip("bx.intervals")

from get_mirna_bs_annotation import site_region

# NM_900001.1 of the synthetic dataset, as gtfToGenePred -genePredExt writes it (0-based starts).
GENEPRED_NM_900001 = "NM_900001.1\tNC_000001.11\t+\t1000\t3400\t1100\t3203\t3\t1000,2000,3000,\t1200,2150,3400,\t0\tSYNA\tcmpl\tcmpl\t0,-1,-1,\n"
# NR_900005.1, the non-coding transcript
GENEPRED_NR_900005 = "NR_900005.1\tNC_000001.11\t+\t12000\t13300\t13300\t13300\t2\t12000,13000,\t12400,13300,\t0\tSYNC\tnone\tnone\t-1,-1,\n"
SITES_HEADER = "FeatureType\tmiRNA\tmRNA\tChrom\tStrand\tG_Start\tG_End\tExon_Starts\tExon_Ends\tLen_NT\tSequence\tSource\tSiteType\n"


def test_refseq_mode_finds_sites_on_the_ensembl_genome(run_script, synthetic_dir, tmp_path):
    """mirna_bs_genomic_coord.py writes sites on Ensembl chromosome "1"; RefSeq isoforms on NC_000001.11 must get them."""
    genepred = tmp_path / "refseq.genepred"
    genepred.write_text(GENEPRED_NM_900001)
    # Bases 11-20 of NM_900001.1 are 1011-1020 of chromosome 1.
    sites = tmp_path / "interactions_gc.txt"
    sites.write_text(SITES_HEADER +
                     "miRNA_binding_site\thsa-miR-1\tNM_900001\t1\t+\t1011\t1020\t1011\t1020\t10\tGCCTGAGTGT\tmiRWalk\t8mer\n")
    output = tmp_path / "mirna_bs_annotation.txt"
    result = run_script("get_mirna_bs_annotation.py", "--chr_ref", synthetic_dir / "chr2acc", "--genepred", genepred,
                        "--isoform_fasta", synthetic_dir / "refseq_rna.fna", "--mirwalk_genomic", sites,
                        "--classification", synthetic_dir / "refseq_sqanti_classification.txt", "--regions", "5UTR",
                        "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == ["NM_900001.1\tmiRWalk\tmiRNA_8mer\t11\t20\t.\t+\t.\tID=hsa-miR-1; Desc=5UTR_8mer"]


@pytest.mark.parametrize("t_start, t_end, cds_start, cds_end, region", [
    (90, 100, "101", "550", "5UTR"),
    (95, 105, "101", "550", "CDS"),
    (540, 560, "101", "550", "CDS"),
    (551, 560, "101", "550", "3UTR"),
    (10, 20, None, None, "noncoding"),
])
def test_site_region(t_start, t_end, cds_start, cds_end, region):
    """Transcript coordinates are already strand-aware; a site that touches the CDS counts as CDS."""
    assert site_region(t_start, t_end, cds_start, cds_end) == region


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        else:
            seqs[name] += line.strip()
    return seqs


@pytest.mark.parametrize("regions, kept", [
    ([], ["3UTR_7mer-A1", "noncoding_6mer"]),
    (["--regions", "5UTR", "CDS", "3UTR", "noncoding"], ["5UTR_8mer", "CDS_7mer-m8", "3UTR_7mer-A1", "noncoding_6mer"]),
])
def test_region_per_isoform(run_script, synthetic_dir, tmp_path, regions, kept):
    """
    Each site gets the region it has on the isoform, from the isoform's CDS (NM_900001.1: 101-550), and only the
    regions asked for are kept: 3UTR and noncoding by default. Lines end in "\\n", not the csv default "\\r\\n".
    """
    genepred = tmp_path / "refseq.genepred"
    genepred.write_text(GENEPRED_NM_900001 + GENEPRED_NR_900005)
    seqs = read_fasta(synthetic_dir / "refseq_rna.fna")
    # (miRNA, transcript, transcript start, end, genomic start, end, site type)
    sites_in = [("hsa-miR-1", "NM_900001.1", 11, 20, 1011, 1020, "8mer"),
                ("hsa-miR-2", "NM_900001.1", 300, 309, 2100, 2109, "7mer-m8"),
                ("hsa-miR-3", "NM_900001.1", 600, 609, 3250, 3259, "7mer-A1"),
                ("hsa-miR-4", "NR_900005.1", 50, 59, 12050, 12059, "6mer")]
    sites = tmp_path / "interactions_gc.txt"
    sites.write_text(SITES_HEADER + "".join(
        f"miRNA_binding_site\t{mirna}\t{tx}\t1\t+\t{g_start}\t{g_end}\t{g_start}\t{g_end}\t10\t{seqs[tx][start - 1:end]}\tmiRWalk\t{kind}\n"
        for mirna, tx, start, end, g_start, g_end, kind in sites_in))
    output = tmp_path / "mirna_bs_annotation.txt"
    result = run_script("get_mirna_bs_annotation.py", "--chr_ref", synthetic_dir / "chr2acc", "--genepred", genepred,
                        "--isoform_fasta", synthetic_dir / "refseq_rna.fna", "--mirwalk_genomic", sites,
                        "--classification", synthetic_dir / "refseq_sqanti_classification.txt", *regions,
                        "--output", output)
    assert result.returncode == 0, result.stderr
    assert b"\r" not in output.read_bytes()
    rows = [line.split("\t") for line in output.read_text().splitlines()]
    assert [row[8].split("Desc=")[1] for row in rows] == kept
    features = {"8mer": "miRNA_8mer", "7mer-m8": "miRNA_7mer", "7mer-A1": "miRNA_7mer", "6mer": "miRNA_6mer"}
    assert [row[2] for row in rows] == [features[desc.split("_")[1]] for desc in kept]
