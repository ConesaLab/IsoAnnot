import io
import os
import sys
from unittest.mock import MagicMock

import pandas as pd
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

import save_biomart_tables
from save_biomart_tables import save_biomart_table

# Nothing listens on port 9: a layer that tried BioMart instead of the saved table would fail
UNREACHABLE_HOST = "http://127.0.0.1:9"

# A saved table as save_biomart_tables.py writes it. SYNA has two gene IDs, SYNB has no GO term,
# and ENSG00000000003 has no gene name.
SAVED_TABLE = (
    "Gene stable ID\tGene name\tGO term accession\tGO term name\tGO domain\n"
    "ENSG00000000001\tSYNA\tGO:0005634\tnucleus\tcellular_component\n"
    "ENSG00000000001\tSYNA\tGO:0003677\tDNA binding\tmolecular_function\n"
    "ENSG00000099999\tSYNA\tGO:0005634\tnucleus\tcellular_component\n"
    "ENSG00000000002\tSYNB\t\t\t\n"
    "ENSG00000000003\t\tGO:0008150\tbiological_process\tbiological_process\n"
)


def write_table(tmp_path, text=SAVED_TABLE):
    table = tmp_path / "biomart_table.tsv"
    table.write_text(text)
    return table


def test_layer_go_reads_saved_table(run_script, synthetic_dir, tmp_path):
    """With --biomart_table, layer_go uses the saved table and never contacts BioMart."""
    output = tmp_path / "layer_go.gtf"
    result = run_script("layer_go.py", "--classification_file", synthetic_dir / "ensembl_sqanti_classification.txt",
                        "--biomart_host", UNREACHABLE_HOST, "--biomart_dataset", "hsapiens_gene_ensembl",
                        "--db", "ensembl", "--biomart_table", write_table(tmp_path), "--output", output)
    assert result.returncode == 0, result.stderr
    assert "Using the saved BioMart table" in result.stderr
    # The second gene ID of SYNA adds no duplicate line, as the unique rows of a BioMart query would
    assert output.read_text().splitlines() == [
        f"{transcript}\tGeneOntology\t{feature}\t.\t.\t.\t.\t.\tID={go}; Name={name}"
        for transcript in ["ENST00000000001", "ENST00000000002"]
        for feature, go, name in [("C", "GO:0005634", "nucleus"), ("F", "GO:0003677", "DNA binding")]
    ]


def test_layer_go_saved_table_mytranscripts_uses_gene_ids(run_script, synthetic_dir, tmp_path):
    """In mytranscripts mode the classification has Ensembl gene IDs, so the saved table is matched on them."""
    classification = tmp_path / "classification.txt"
    classification.write_text((synthetic_dir / "ensembl_sqanti_classification.txt").read_text().replace("\tSYNA\t", "\tENSG00000099999\t"))
    output = tmp_path / "layer_go.gtf"
    result = run_script("layer_go.py", "--classification_file", classification,
                        "--biomart_host", UNREACHABLE_HOST, "--biomart_dataset", "hsapiens_gene_ensembl",
                        "--db", "mytranscripts", "--biomart_table", write_table(tmp_path), "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [
        "ENST00000000001\tGeneOntology\tC\t.\t.\t.\t.\t.\tID=GO:0005634; Name=nucleus",
        "ENST00000000002\tGeneOntology\tC\t.\t.\t.\t.\t.\tID=GO:0005634; Name=nucleus",
    ]


def test_layer_reactome_reads_saved_table(run_script, synthetic_dir, tmp_path):
    """With --biomart_table, layer_reactome takes gene ID and name pairs from the saved table."""
    reactome = tmp_path / "Ensembl2Reactome.txt"
    reactome.write_text(
        "ENSG00000000001\tR-HSA-1\thttps://reactome.org/R-HSA-1\tPathway one\tTAS\tHomo sapiens\n"
        "ENSG00000000003\tR-HSA-3\thttps://reactome.org/R-HSA-3\tPathway three\tTAS\tHomo sapiens\n"
        "ENSG00000000001\tR-MMU-1\thttps://reactome.org/R-MMU-1\tMouse pathway\tTAS\tMus musculus\n"
    )
    output = tmp_path / "layer_reactome.gtf"
    result = run_script("layer_reactome.py", "--reactome_file", reactome,
                        "--classification_file", synthetic_dir / "ensembl_sqanti_classification.txt",
                        "--biomart_host", UNREACHABLE_HOST, "--biomart_dataset", "hsapiens_gene_ensembl",
                        "--species", "Homo sapiens", "--db", "ensembl", "--biomart_table", write_table(tmp_path),
                        "--output", output)
    assert result.returncode == 0, result.stderr
    # ENSG00000000003 has no gene name, so its pathway is not attached, as with a BioMart query
    assert output.read_text().splitlines() == [
        "ENST00000000001\tREACTOME\tPATHWAY\t.\t.\t.\t.\t.\tID=R-HSA-1; Name=Pathway one",
        "ENST00000000002\tREACTOME\tPATHWAY\t.\t.\t.\t.\t.\tID=R-HSA-1; Name=Pathway one",
    ]


def test_layer_go_rejects_a_table_that_is_not_saved_by_isoannot(run_script, synthetic_dir, tmp_path):
    """A table without the BioMart columns stops the layer with a message that names the script to use."""
    table = write_table(tmp_path, "gene\tgo\nSYNA\tGO:0005634\n")
    result = run_script("layer_go.py", "--classification_file", synthetic_dir / "ensembl_sqanti_classification.txt",
                        "--biomart_host", UNREACHABLE_HOST, "--biomart_dataset", "hsapiens_gene_ensembl",
                        "--db", "ensembl", "--biomart_table", table, "--output", tmp_path / "layer_go.gtf")
    assert result.returncode != 0
    assert "is not a saved BioMart table: columns Gene stable ID, Gene name, GO term accession, GO term name, GO domain are missing" in result.stderr


def test_save_biomart_table_writes_what_the_layers_read(tmp_path):
    """save_biomart_tables.py writes BioMart's answer with its own column names, in the format read_biomart_table expects."""
    answer = pd.DataFrame([["ENSG00000000001", "SYNA", "GO:0005634", "nucleus", "cellular_component"],
                           ["ENSG00000000002", "SYNB", None, None, None]],
                          columns=["Gene stable ID", "Gene name", "GO term accession", "GO term name", "GO domain"])
    dataset = MagicMock()
    dataset.query.return_value = answer
    output = tmp_path / "biomart_table.tsv"

    assert save_biomart_table(dataset, str(output)) == 2
    assert dataset.query.call_args.kwargs["attributes"] == ["ensembl_gene_id", "external_gene_name", "go_id", "name_1006", "namespace_1003"]
    assert output.read_text().splitlines() == [
        "Gene stable ID\tGene name\tGO term accession\tGO term name\tGO domain",
        "ENSG00000000001\tSYNA\tGO:0005634\tnucleus\tcellular_component",
        "ENSG00000000002\tSYNB\t\t\t",
    ]
    assert not (tmp_path / "biomart_table.tsv.tmp").exists()


def test_https_dataset_reads_a_complete_answer(monkeypatch):
    """Archive hosts only answer over https, which pybiomart 0.2.0 can't use, so the table is saved over https."""
    calls = []

    def fake_urlopen(url, timeout):
        calls.append(url)
        return io.BytesIO(b"Gene stable ID\tGene name\nENSG00000000001\tSYNA\n[success]\n")
    monkeypatch.setattr(save_biomart_tables.urllib.request, "urlopen", fake_urlopen)

    table = save_biomart_tables.https_dataset("sep2025.archive.ensembl.org", "hsapiens_gene_ensembl").query(["ensembl_gene_id", "external_gene_name"])
    assert calls[0].startswith("https://sep2025.archive.ensembl.org/biomart/martservice?query=")
    assert 'completionStamp%3D%221%22' in calls[0]
    assert table.values.tolist() == [["ENSG00000000001", "SYNA"]]


def test_https_dataset_rejects_an_incomplete_answer(monkeypatch):
    """BioMart can stop half-way without an error; an answer without the [success] line is not saved."""
    monkeypatch.setattr(save_biomart_tables.urllib.request, "urlopen",
                        lambda url, timeout: io.BytesIO(b"Gene stable ID\tGene name\nENSG00000000001\tSY"))

    with pytest.raises(ValueError, match=r"is incomplete \(no \[success\] line at the end\)"):
        save_biomart_tables.https_dataset("sep2025.archive.ensembl.org", "hsapiens_gene_ensembl").query(["ensembl_gene_id", "external_gene_name"])
