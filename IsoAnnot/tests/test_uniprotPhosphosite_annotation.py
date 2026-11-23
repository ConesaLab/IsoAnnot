import os
import sys

import pytest

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


def test_matched_long_read_orf_does_not_reuse_the_position_of_the_previous_motif(run_script, synthetic_dir, tmp_path):
    """
    The second motif of the matched protein lies in an intron of PB.1.1. It used to be written at the
    position of the motif before it (serine 6), because the residue there is also a serine.
    """
    lines = run_long_read_isoform(run_script, synthetic_dir, tmp_path,
                                  [motif_line("MOD_RES", "Phosphoserine", 1116, 1118, "S", "NP_900001"),
                                   motif_line("MOD_RES", "N-acetylserine", 1500, 1502, "S", "NP_900001")],
                                  ["PB.1.1\tP0SYA1\tswissprot\tsequence\n", "PB.1.1\tNP_900001\trefseq\tsequence\n"])
    assert lines == [PHOSPHOSERINE_6]


def test_matched_long_read_orf_with_versioned_refseq_id(run_script, synthetic_dir, tmp_path):
    """
    protein_assoc has the RefSeq IDs with their version and the motifs don't (NP_900001.1 and NP_900001):
    the ORF still gets the motifs of its protein.
    """
    lines = run_long_read_isoform(run_script, synthetic_dir, tmp_path,
                                  [motif_line("MOD_RES", "Phosphoserine", 1116, 1118, "S", "NP_900001")],
                                  ["PB.1.1\tP0SYA1\tswissprot\tsequence\n", "PB.1.1\tNP_900001.1\trefseq\tsequence\n"])
    assert lines == [PHOSPHOSERINE_6]


def run_refseq_motifs(run_script, synthetic_dir, tmp_path, motif_lines, *options):
    """Runs the script in refseq mode on NM_900001.1 with the given motifs, and returns the output lines and the log."""
    genepred = tmp_path / "refseq.genepred"
    genepred.write_text(GENEPRED_NM_900001)
    motifs = tmp_path / "uniprot_gcord_domainGenomic.txt"
    motifs.write_text("".join(motif_lines))
    output = tmp_path / "uniprot_Phosphosite_info.txt"
    result = run_script("uniprotPhosphosite_annotation.py", "--orf_fasta", synthetic_dir / "refseq_protein.faa",
                        "--classification_file", synthetic_dir / "refseq_sqanti_classification.txt",
                        "--genepred_file", genepred, "--uniprotmotif_file", motifs,
                        "--chr_ref", synthetic_dir / "chr2acc", "--keep_version", "True", "--db", "refseq",
                        "--biomart_host", "www.ensembl.org", "--output", output, *options)
    assert result.returncode == 0, result.stderr
    return output.read_text().splitlines(), result.stderr


def test_note_in_brackets_before_the_residue_is_not_part_of_the_id(run_script, synthetic_dir, tmp_path):
    """tappAS rejects IDs with spaces: "(Microbial infection) Phosphoserine" was written whole as the ID."""
    lines, log = run_refseq_motifs(run_script, synthetic_dir, tmp_path,
                                   [motif_line("MOD_RES", "(Microbial infection) Phosphoserine; by host", 1116, 1118, "S", "NP_900001")])
    assert lines == ["NM_900001.1\tUniProtKB/Swiss-Prot_Phosphosite\tPTM\t6\t6\t.\t.\t.\t"
                     "ID=Phosphoserine; Name=Phosphoserine; Desc=(Microbial Infection) Phosphoserine_SwissProt"]
    assert "ID with spaces" not in log


def test_any_other_id_with_spaces_is_joined_and_reported_once(run_script, synthetic_dir, tmp_path):
    """A wording that no rule covers still gives an ID tappAS accepts, and one warning with the count."""
    lines, log = run_refseq_motifs(run_script, synthetic_dir, tmp_path,
                                   [motif_line("MOD_RES", "Phosphoserine or phosphothreonine; alternate", 1116, 1118, "S", "NP_900001")])
    assert lines == ["NM_900001.1\tUniProtKB/Swiss-Prot_Phosphosite\tPTM\t6\t6\t.\t.\t.\t"
                     "ID=Phosphoserine_Or_Phosphothreonine; Name=Phosphoserine_Or_Phosphothreonine; "
                     "Desc=Phosphoserine Or Phosphothreonine_SwissProt"]
    assert log.count("ID with spaces") == 1 and "1 UniProt features had an ID with spaces" in log


def psp_line(description, start, end, sequence, protein_id):
    """One PhosphoSitePlus line of uniprot_gcord_domainGenomic.txt, as uniprotPhosphosite_genomicCoordinates.py writes it."""
    return "\t".join(["MOD_RES", description, "1", "+", str(start), str(end), str(start), str(end), "1", "3", sequence,
                      "PhosphositePlus", "P0SYA1", protein_id, f"{sequence}6-p"]) + "\n"


# Serine 6 of NP_900001.1 (codon 1116-1118) with three modifications: a phosphorylation in both sources,
# an acetylation only in Swiss-Prot and a methylation whose degree differs between the sources
BOTH_SOURCES = [motif_line("MOD_RES", "Phosphoserine; by PKA", 1116, 1118, "S", "NP_900001"),
                psp_line("phosphoserine", 1116, 1118, "S", "NP_900001"),
                motif_line("MOD_RES", "N-acetylserine", 1116, 1118, "S", "NP_900001"),
                motif_line("MOD_RES", "N6,N6-dimethylserine", 1116, 1118, "S", "NP_900001"),
                psp_line("methylserine", 1116, 1118, "S", "NP_900001")]
ROW = "NM_900001.1\tUniProtKB/Swiss-Prot_Phosphosite\tPTM\t6\t6\t.\t.\t.\t"
SWISSPROT_ROWS = [ROW + "ID=Phosphoserine; Name=Phosphoserine; Desc=Phosphoserine_SwissProt",
                  ROW + "ID=Acetylation; Name=Acetylation; Desc=N-Acetylserine_SwissProt",
                  ROW + "ID=Methylation; Name=Methylation; Desc=N6,N6-Dimethylserine_SwissProt"]
METHYLSERINE = ROW + "ID=Methylation; Name=Methylation; Desc=Methylserine_PhosphositePlus"


def test_ptm_reported_by_both_sources_is_written_once_as_swissprot(run_script, synthetic_dir, tmp_path):
    """
    The phosphoserine of PhosphoSitePlus is left out, because Swiss-Prot has it. Its methylserine stays:
    the Swiss-Prot row at that residue is a dimethylation, which is another modification.
    """
    lines, log = run_refseq_motifs(run_script, synthetic_dir, tmp_path, BOTH_SOURCES)
    assert sorted(lines) == sorted(SWISSPROT_ROWS + [METHYLSERINE])
    assert "1 PhosphoSitePlus modifications were left out" in log


def test_keep_both_sources_writes_the_ptm_once_for_each_source(run_script, synthetic_dir, tmp_path):
    lines, log = run_refseq_motifs(run_script, synthetic_dir, tmp_path, BOTH_SOURCES, "--keep_both_sources")
    assert sorted(lines) == sorted(SWISSPROT_ROWS + [METHYLSERINE, ROW + "ID=Phosphoserine; Name=Phosphoserine; Desc=Phosphoserine_PhosphositePlus"])


GENEPRED_ENST_1 = GENEPRED_NM_900001.replace("NM_900001.1\tNC_000001.11", "ENST00000000001\t1")


def run_ensembl(run_script, synthetic_dir, tmp_path, genepred_text, keep_version="False", classification=None):
    """Runs the script in Ensembl mode with a phosphoserine on serine 6 of ENSP00000000001."""
    genepred = tmp_path / "ensembl.genepred"
    genepred.write_text(genepred_text)
    motifs = tmp_path / "uniprot_gcord_domainGenomic.txt"
    motifs.write_text(motif_line("MOD_RES", "Phosphoserine", 1116, 1118, "S", "ENSP00000000001"))
    output = tmp_path / "uniprot_Phosphosite_info.txt"
    result = run_script("uniprotPhosphosite_annotation.py", "--orf_fasta", synthetic_dir / "ensembl_pep.fa",
                        "--classification_file", classification or synthetic_dir / "ensembl_sqanti_classification.txt",
                        "--genepred_file", genepred, "--uniprotmotif_file", motifs,
                        "--chr_ref", synthetic_dir / "chr2acc", "--keep_version", keep_version, "--db", "ensembl",
                        "--biomart_host", "e115.ensembl.org", "--output", output)
    return result, output


@pytest.mark.parametrize("keep_version", ["False", "True"])
def test_ensembl_protein_is_found_whatever_transcript_versioned_says(run_script, synthetic_dir, tmp_path, keep_version):
    """
    The Ensembl protein FASTA has versions (ENSP00000000001.2) and the classification does not. With
    transcript_versioned: True in an Ensembl config no transcript found its protein, and the layer was empty.
    """
    result, output = run_ensembl(run_script, synthetic_dir, tmp_path, GENEPRED_ENST_1, keep_version)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [PHOSPHOSERINE_6.replace("PB.1.1", "ENST00000000001")]
    assert "the UniProt layer will be empty" not in result.stderr


def test_warns_when_no_coding_transcript_has_its_protein(run_script, synthetic_dir, tmp_path):
    """A classification and a protein FASTA of different annotations give an empty layer: it is said once."""
    classification = tmp_path / "classification.txt"
    classification.write_text((synthetic_dir / "ensembl_sqanti_classification.txt").read_text().replace("\tENSP0000", "\tOTHERP0000"))
    result, output = run_ensembl(run_script, synthetic_dir, tmp_path, GENEPRED_ENST_1, classification=classification)
    assert result.returncode == 0, result.stderr
    assert output.read_text() == ""
    assert result.stderr.count("coding transcripts of") == 1
    assert "has its protein in the protein FASTA, so the UniProt layer will be empty" in result.stderr


def test_stops_when_the_table_does_not_know_the_chromosomes_of_the_transcripts(run_script, synthetic_dir, tmp_path):
    """Transcripts on chromosomes that chr2acc can't translate got no motif, without an error."""
    result, output = run_ensembl(run_script, synthetic_dir, tmp_path, GENEPRED_ENST_1.replace("\t1\t+", "\tNW_999999.1\t+"))
    assert result.returncode != 0
    assert "Critical Error in uniprotPhosphosite_annotation.py: 1/1 chromosome lookups failed" in result.stderr
