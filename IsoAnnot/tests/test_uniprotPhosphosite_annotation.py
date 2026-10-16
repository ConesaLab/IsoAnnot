import os
import sys

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

# NM_900001.1 of the synthetic dataset, as gtfToGenePred -genePredExt writes it (0-based starts).
GENEPRED_NM_900001 = "NM_900001.1\tNC_000001.11\t+\t1000\t3400\t1100\t3203\t3\t1000,2000,3000,\t1200,2150,3400,\t0\tSYNA\tcmpl\tcmpl\t0,-1,-1,\n"


def test_refseq_mode_projects_motifs_from_the_ensembl_genome(run_script, synthetic_dir, tmp_path):
    """The motif file is on Ensembl chromosome "1"; RefSeq transcripts on NC_000001.11 must still get the motif."""
    genepred = tmp_path / "refseq.genepred"
    genepred.write_text(GENEPRED_NM_900001)
    # Serine 6 of NP_900001.1 (MDQLNS...) is the codon 1116-1118 of chromosome 1, in the first CDS exon.
    motifs = tmp_path / "uniprot_gcord_domainGenomic.txt"
    motifs.write_text("\t".join(["MOD_RES", "Phosphoserine", "1", "+", ".", ".", "1116", "1118", ".", ".", "S",
                                 "UniProtKB", "P0SYA1", "NP_900001"]) + "\n")
    output = tmp_path / "uniprot_Phosphosite_info.txt"
    result = run_script("uniprotPhosphosite_annotation.py", "--orf_fasta", synthetic_dir / "refseq_protein.faa",
                        "--classification_file", synthetic_dir / "refseq_sqanti_classification.txt",
                        "--genepred_file", genepred, "--uniprotmotif_file", motifs,
                        "--chr_ref", synthetic_dir / "chr2acc", "--keep_version", "True", "--db", "refseq",
                        "--biomart_host", "https://www.ensembl.org", "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [
        "NM_900001.1\tUniProtKB/Swiss-Prot_Phosphosite\tPTM\t6\t6\t.\t.\t.\t"
        "ID=Phosphoserine; Name=Phosphoserine; Desc=Phosphoserine_UniProtKB"]
