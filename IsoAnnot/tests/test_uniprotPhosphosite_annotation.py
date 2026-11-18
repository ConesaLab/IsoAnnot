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


def motif_line(feature, description, start, end, sequence, protein_id, strand="+"):
    """One line of uniprot_gcord_domainGenomic.txt, on Ensembl chromosome 1."""
    return "\t".join([feature, description, "1", strand, ".", ".", str(start), str(end), ".", ".", sequence,
                      "SwissProt", "P0SYA1", protein_id]) + "\n"


def run_long_read_isoform(run_script, synthetic_dir, tmp_path, motif_lines, assoc_lines):
    """
    Runs the script in mytranscripts mode on PB.1.1, a long-read isoform with the exons and the ORF of
    NM_900001.1 (protein MDQLNSAS...), and returns the output lines.
    """
    genepred = tmp_path / "mytranscripts.genepred"
    genepred.write_text(GENEPRED_NM_900001.replace("NM_900001.1\tNC_000001.11", "PB.1.1\t1"))
    classification = tmp_path / "sqanti_classification.txt"
    classification.write_text("isoform\tcoding\tCDS_start\tCDS_end\nPB.1.1\tcoding\t101\t550\n")
    protein = [line for line in (synthetic_dir / "refseq_protein.faa").read_text().split(">") if line.startswith("NP_900001.1")][0]
    orfs = tmp_path / "final_sqanti_corrected.faa"
    orfs.write_text(">PB.1.1 synthetic ORF\n" + protein.split("\n", 1)[1])
    motifs = tmp_path / "uniprot_gcord_domainGenomic.txt"
    motifs.write_text("".join(motif_lines))
    assoc = tmp_path / "protein_assoc_data.txt"
    assoc.write_text("".join(assoc_lines))
    output = tmp_path / "uniprot_Phosphosite_info.txt"
    result = run_script("uniprotPhosphosite_annotation.py", "--orf_fasta", orfs, "--classification_file", classification,
                        "--genepred_file", genepred, "--uniprotmotif_file", motifs, "--chr_ref", synthetic_dir / "chr2acc",
                        "--protein_association", assoc, "--keep_version", "True", "--db", "mytranscripts",
                        "--biomart_host", "www.ensembl.org", "--output", output)
    assert result.returncode == 0, result.stderr
    return output.read_text().splitlines()


PHOSPHOSERINE_6 = ("PB.1.1\tUniProtKB/Swiss-Prot_Phosphosite\tPTM\t6\t6\t.\t.\t.\t"
                   "ID=Phosphoserine; Name=Phosphoserine; Desc=Phosphoserine_SwissProt")


def test_reference_mode_only_gives_a_transcript_the_motifs_of_its_protein(run_script, synthetic_dir, tmp_path):
    """Motifs are looked up by protein ID: the serine of another protein at the same codon is not written."""
    genepred = tmp_path / "refseq.genepred"
    genepred.write_text(GENEPRED_NM_900001)
    motifs = tmp_path / "uniprot_gcord_domainGenomic.txt"
    motifs.write_text(motif_line("MOD_RES", "N-acetylserine", 1116, 1118, "S", "NP_900009")
                      + motif_line("MOD_RES", "Phosphoserine", 1116, 1118, "S", "NP_900001"))
    output = tmp_path / "uniprot_Phosphosite_info.txt"
    result = run_script("uniprotPhosphosite_annotation.py", "--orf_fasta", synthetic_dir / "refseq_protein.faa",
                        "--classification_file", synthetic_dir / "refseq_sqanti_classification.txt",
                        "--genepred_file", genepred, "--uniprotmotif_file", motifs,
                        "--chr_ref", synthetic_dir / "chr2acc", "--keep_version", "True", "--db", "refseq",
                        "--biomart_host", "www.ensembl.org", "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [PHOSPHOSERINE_6.replace("PB.1.1", "NM_900001.1")]


def test_novel_long_read_orf_gets_the_motifs_it_overlaps(run_script, synthetic_dir, tmp_path):
    """
    A novel ORF has no protein to look up, so its motifs are found by position: the serine inside its
    CDS is written, the one 5 Mb away (another genomic bin) and the one on the other strand are not.
    """
    lines = run_long_read_isoform(run_script, synthetic_dir, tmp_path,
                                  [motif_line("MOD_RES", "Phosphoserine", 1116, 1118, "S", "ENSP_OTHER"),
                                   motif_line("MOD_RES", "Phosphoserine", 1116, 1118, "S", "ENSP_MINUS", strand="-"),
                                   motif_line("MOD_RES", "Phosphoserine", 5001116, 5001118, "S", "ENSP_FAR")],
                                  ["PB.1.1\tProtein_PB.1.1\tnovel\tsequence\n"])
    assert lines == [PHOSPHOSERINE_6]
