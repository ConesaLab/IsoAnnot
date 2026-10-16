import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

pytest.importorskip("bx.intervals")

# NM_900001.1 of the synthetic dataset, as gtfToGenePred -genePredExt writes it (0-based starts).
GENEPRED_NM_900001 = "NM_900001.1\tNC_000001.11\t+\t1000\t3400\t1100\t3203\t3\t1000,2000,3000,\t1200,2150,3400,\t0\tSYNA\tcmpl\tcmpl\t0,-1,-1,\n"


def test_refseq_mode_finds_sites_on_the_ensembl_genome(run_script, synthetic_dir, tmp_path):
    """mirna_bs_genomic_coord.py writes sites on Ensembl chromosome "1"; RefSeq isoforms on NC_000001.11 must get them."""
    genepred = tmp_path / "refseq.genepred"
    genepred.write_text(GENEPRED_NM_900001)
    # Bases 11-20 of NM_900001.1 are 1011-1020 of chromosome 1.
    sites = tmp_path / "interactions_gc.txt"
    sites.write_text("FeatureType\tmiRNA\tmRNA\tChrom\tStrand\tG_Start\tG_End\tExon_Starts\tExon_Ends\tLen_NT\tSequence\tSource\n"
                     "miRNA_binding_site\thsa-miR-1\tNM_900001\t1\t+\t1011\t1020\t1011\t1020\t10\tGCCTGAGTGT\tmiRWalk\n")
    output = tmp_path / "mirna_bs_annotation.txt"
    result = run_script("get_mirna_bs_annotation.py", "--chr_ref", synthetic_dir / "chr2acc", "--genepred", genepred,
                        "--isoform_fasta", synthetic_dir / "refseq_rna.fna", "--mirwalk_genomic", sites,
                        "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == ["NM_900001.1\tmiRWalk\tmiRNA_binding_site\t11\t20\t.\t+\t.\tID=hsa-miR-1"]
