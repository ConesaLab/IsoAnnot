import json
import logging
import os
import sys

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from IsoAnnot import annotation_id_forms, read_gtf_protein_ids


def test_annotation_id_forms_ensembl_drops_the_version(caplog):
    """
    The Ensembl GTF has the IDs without version, so they are stored without it. A protein that is not in
    the annotation (human proteins on scaffolds, which the .chr GTF leaves out) follows the same form.
    """
    with caplog.at_level(logging.INFO):
        forms = annotation_id_forms(["ENSP00000397289.1", "ENSP00000216124.10", "ENSP00000487140.1"],
                                    {"ENSP00000397289", "ENSP00000216124"}, "ensembl")
    assert forms == {"ENSP00000397289.1": "ENSP00000397289", "ENSP00000216124.10": "ENSP00000216124",
                     "ENSP00000487140.1": "ENSP00000487140"}
    assert [record.levelname for record in caplog.records] == ["INFO"]
    assert "0 protein IDs are written as in the annotation, 2 without their final number" in caplog.text
    assert "1 are not in the annotation" in caplog.text


def test_annotation_id_forms_refseq_keeps_the_version(caplog):
    """The RefSeq GTF has versioned IDs, so the version stays, also for a protein that is not in the GTF."""
    with caplog.at_level(logging.INFO):
        forms = annotation_id_forms(["NP_000005.3", "XP_047292308.1", "NP_000238.1"],
                                    {"NP_000005.3", "XP_047292308.1"}, "refseq")
    assert forms == {"NP_000005.3": "NP_000005.3", "XP_047292308.1": "XP_047292308.1", "NP_000238.1": "NP_000238.1"}
    assert [record.levelname for record in caplog.records] == ["INFO"]


def test_annotation_id_forms_isoform_numbers():
    """In Ensembl Plants the number is the isoform and the GTF has it, so two isoforms stay two IDs."""
    forms = annotation_id_forms(["AT1G01010.1", "AT1G01010.2", "FBpp0070000", "PB.1.1"],
                                {"AT1G01010.1", "AT1G01010.2", "FBpp0070000"}, "ensembl")
    assert forms == {"AT1G01010.1": "AT1G01010.1", "AT1G01010.2": "AT1G01010.2", "FBpp0070000": "FBpp0070000",
                     "PB.1.1": "PB.1.1"}


def test_annotation_id_forms_mixed_forms(caplog):
    """With both forms in one database there is no form to follow: the other IDs stay as in the FASTA."""
    forms = annotation_id_forms(["ENSP00000000001.2", "NP_900001.1", "NEW_1.4"], {"ENSP00000000001", "NP_900001.1"}, "custom")
    assert forms == {"ENSP00000000001.2": "ENSP00000000001", "NP_900001.1": "NP_900001.1", "NEW_1.4": "NEW_1.4"}
    assert "do not share one form of protein ID (1 IDs are the same, 1 differ in the final number)" in caplog.text
    assert "The 1 protein IDs that are not in the annotation are written as in the protein FASTA" in caplog.text


def test_annotation_id_forms_no_id_in_common(caplog):
    """An annotation without the protein IDs of the FASTA changes nothing and says so."""
    forms = annotation_id_forms(["PROT_1.1", "PROT_2.1"], set(), "custom")
    assert forms == {"PROT_1.1": "PROT_1.1", "PROT_2.1": "PROT_2.1"}
    assert "(0 IDs are the same, 0 differ in the final number)" in caplog.text


def test_annotation_id_forms_never_merges_two_ids(caplog):
    """Two IDs outside the annotation that only differ in the number are two proteins: they keep it."""
    forms = annotation_id_forms(["ENSP00000000001.2", "GENE9.1", "GENE9.2", "ENSP00000000009.1"], {"ENSP00000000001"}, "ensembl")
    assert forms == {"ENSP00000000001.2": "ENSP00000000001", "GENE9.1": "GENE9.1", "GENE9.2": "GENE9.2",
                     "ENSP00000000009.1": "ENSP00000000009"}
    assert "2 protein IDs that are not in the annotation keep their final number" in caplog.text


def test_read_gtf_protein_ids(synthetic_dir):
    """The protein IDs are read from the CDS lines as each annotation writes them."""
    assert read_gtf_protein_ids(str(synthetic_dir / "ensembl.gtf")) == {
        "ENSP00000000001", "ENSP00000000002", "ENSP00000000003", "ENSP00000000004", "ENSP00000000006"}
    assert read_gtf_protein_ids(str(synthetic_dir / "refseq.gtf")) == {
        "NP_900001.1", "NP_900002.1", "NP_900003.1", "NP_900004.1", "NP_900006.2"}


def test_parse_protein_databases_id_forms(synthetic_dir, run_script, tmp_path):
    """
    parsed_databases.json has every ID in the form of its annotation. A scaffold protein with the sequence
    of an annotated one is listed without version, like the other Ensembl IDs, and Ensembl Plants IDs keep
    their isoform number.
    """
    pep = tmp_path / "pep.fa"
    pep.write_text(">ENSP00000000001.2 pep chromosome:SYNTH1:1:1101:3200:1 gene:ENSG00000000001.3\nMAAAK\n"
                   ">ENSP00000000099.1 pep scaffold:SYNTH1:HSCHR1_1_CTG3:1:300:1 gene:ENSG00000000099.1\nMAAAK\n")
    plants_pep = tmp_path / "plants_pep.fa"
    plants_pep.write_text(">AT1G01010.1 pep chromosome:TAIR10:1:3631:5899:1 gene:AT1G01010\nMCCCK\n"
                          ">AT1G01010.2 pep chromosome:TAIR10:1:3631:5899:1 gene:AT1G01010\nMGGGK\n")
    plants_gtf = tmp_path / "plants.gtf"
    plants_gtf.write_text("".join(
        f'1\taraport11\tCDS\t3760\t3913\t.\t+\t0\tgene_id "AT1G01010"; transcript_id "{protein_id}"; protein_id "{protein_id}";\n'
        for protein_id in ("AT1G01010.1", "AT1G01010.2")))

    output = tmp_path / "parsed_databases.json"
    result = run_script("parse_protein_databases.py", "--uniprot_fasta", synthetic_dir / "uniprot.fasta",
                        "--ensembl_fasta", pep, "--ensembl_gtf", synthetic_dir / "ensembl.gtf",
                        "--refseq_fasta", synthetic_dir / "refseq_protein.faa", "--refseq_gtf", synthetic_dir / "refseq.gtf",
                        "--output", output)
    assert result.returncode == 0, result.stderr
    parsed = json.loads(output.read_text())
    assert parsed["MAAAK"] == {"ensembl": ["ENSP00000000001", "ENSP00000000099"]}
    assert sorted(protein_id for dbs in parsed.values() for protein_id in dbs.get("refseq", [])) == [
        "NP_900001.1", "NP_900002.1", "NP_900003.1", "NP_900004.1", "NP_900006.2"]
    assert "WARNING" not in result.stderr

    result = run_script("parse_protein_databases.py", "--uniprot_fasta", synthetic_dir / "uniprot.fasta",
                        "--ensembl_fasta", plants_pep, "--ensembl_gtf", plants_gtf, "--output", output)
    assert result.returncode == 0, result.stderr
    parsed = json.loads(output.read_text())
    assert parsed["MCCCK"] == {"ensembl": ["AT1G01010.1"]}
    assert parsed["MGGGK"] == {"ensembl": ["AT1G01010.2"]}
