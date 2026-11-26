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
    """layer_go reads the gene names and GO terms from the saved table."""
    output = tmp_path / "layer_go.gtf"
    result = run_script("layer_go.py", "--classification_file", synthetic_dir / "ensembl_sqanti_classification.txt",
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
                        "--db", "mytranscripts", "--biomart_table", write_table(tmp_path), "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [
        "ENST00000000001\tGeneOntology\tC\t.\t.\t.\t.\t.\tID=GO:0005634; Name=nucleus",
        "ENST00000000002\tGeneOntology\tC\t.\t.\t.\t.\t.\tID=GO:0005634; Name=nucleus",
    ]


def test_layer_reactome_reads_saved_table(run_script, synthetic_dir, tmp_path):
    """layer_reactome takes gene ID and name pairs from the saved table."""
    reactome = tmp_path / "Ensembl2Reactome.txt"
    reactome.write_text(
        "ENSG00000000001\tR-HSA-1\thttps://reactome.org/R-HSA-1\tPathway one\tTAS\tHomo sapiens\n"
        "ENSG00000000003\tR-HSA-3\thttps://reactome.org/R-HSA-3\tPathway three\tTAS\tHomo sapiens\n"
        "ENSG00000000001\tR-MMU-1\thttps://reactome.org/R-MMU-1\tMouse pathway\tTAS\tMus musculus\n"
    )
    output = tmp_path / "layer_reactome.gtf"
    result = run_script("layer_reactome.py", "--reactome_file", reactome,
                        "--classification_file", synthetic_dir / "ensembl_sqanti_classification.txt",
                        "--species", "Homo sapiens", "--db", "ensembl", "--biomart_table", write_table(tmp_path),
                        "--output", output)
    assert result.returncode == 0, result.stderr
    # ENSG00000000003 has no gene name, so its pathway is not attached
    assert output.read_text().splitlines() == [
        "ENST00000000001\tREACTOME\tPATHWAY\t.\t.\t.\t.\t.\tID=R-HSA-1; Name=Pathway one",
        "ENST00000000002\tREACTOME\tPATHWAY\t.\t.\t.\t.\t.\tID=R-HSA-1; Name=Pathway one",
    ]


def test_layer_go_rejects_a_table_that_is_not_saved_by_isoannot(run_script, synthetic_dir, tmp_path):
    """A table without the BioMart columns stops the layer with a message that names the script to use."""
    table = write_table(tmp_path, "gene\tgo\nSYNA\tGO:0005634\n")
    result = run_script("layer_go.py", "--classification_file", synthetic_dir / "ensembl_sqanti_classification.txt",
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
    """Archive hosts only answer over https, so the table is asked for over https."""
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


@pytest.mark.parametrize("layer, options", [
    ("layer_go.py", []),
    ("layer_reactome.py", ["--reactome_file", "Ensembl2Reactome.txt", "--species", "Homo sapiens"]),
])
def test_layers_need_the_saved_table(run_script, synthetic_dir, tmp_path, layer, options):
    """The layers no longer ask BioMart themselves: without a table they stop at the arguments."""
    result = run_script(layer, "--classification_file", synthetic_dir / "ensembl_sqanti_classification.txt",
                        "--db", "ensembl", "--output", tmp_path / "layer.gtf", *options)
    assert result.returncode != 0
    assert "the following arguments are required: --biomart_table" in result.stderr


@pytest.mark.parametrize("host", ["e115.ensembl.org", "http://e115.ensembl.org", "https://e115.ensembl.org", "https://e115.ensembl.org/", " e115.ensembl.org"])
def test_biomart_url_accepts_the_host_in_any_form(host):
    """A host with https:// used to fail in the layers, and one with http:// is refused by the archives."""
    assert save_biomart_tables.biomart_url(host) == "https://e115.ensembl.org/biomart/martservice"


@pytest.mark.parametrize("registry, release", [
    ('<MartRegistry>\n  <MartURLLocation database="ensembl_mart_115" default="1" displayName="Ensembl Genes 115" host="sep2025.archive.ensembl.org" />\n</MartRegistry>', 115),
    ('<MartRegistry>\n  <MartURLLocation database="plants_mart_62" default="1" displayName="Ensembl Plants Genes 62" />\n</MartRegistry>', 62),
    ("<html>Service unavailable</html>", None),
])
def test_biomart_release_reads_the_gene_database_name(monkeypatch, registry, release):
    calls = []

    def fake_urlopen(url, timeout):
        calls.append(url)
        return io.BytesIO(registry.encode())
    monkeypatch.setattr(save_biomart_tables.urllib.request, "urlopen", fake_urlopen)

    assert save_biomart_tables.biomart_release("e115.ensembl.org") == release
    assert calls == ["https://e115.ensembl.org/biomart/martservice?type=registry"]


def test_biomart_release_is_unknown_when_the_host_does_not_answer(monkeypatch):
    """The release check is only a warning, so a host that does not answer it does not stop anything."""
    def fake_urlopen(url, timeout):
        raise OSError("Name or service not known")
    monkeypatch.setattr(save_biomart_tables.urllib.request, "urlopen", fake_urlopen)

    assert save_biomart_tables.biomart_release("e115.ensembl.org") is None


@pytest.mark.parametrize("url, release", [
    ("https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.chr.gtf.gz", 115),
    ("https://ftp.ensemblgenomes.org/pub/release-62/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gtf.gz", 62),
    ("/home/user/my_annotation.gtf", None),
    (None, None),
])
def test_release_of_url(url, release):
    assert save_biomart_tables.release_of_url(url) == release


@pytest.mark.parametrize("biomart, files, warns", [
    (116, 115, True),
    (115, 115, False),
    (None, 115, False),
    (116, None, False),
])
def test_release_mismatch_warns_only_when_both_releases_are_known_and_differ(biomart, files, warns):
    warning = save_biomart_tables.release_mismatch(biomart, files)
    assert (warning is not None) == warns
    if warns:
        assert warning.startswith("BioMart serves Ensembl release 116, but the files of the config are release 115.")
        assert "e115.ensembl.org" in warning


def run_main(monkeypatch, tmp_path, registry_release, *options):
    """Runs save_biomart_tables.main() with a BioMart that serves registry_release and answers one gene."""
    def fake_urlopen(url, timeout):
        if url.endswith("?type=registry"):
            return io.BytesIO(f'<MartURLLocation database="ensembl_mart_{registry_release}" />'.encode())
        return io.BytesIO(b"Gene stable ID\tGene name\tGO term accession\tGO term name\tGO domain\n"
                          b"ENSG00000000001\tSYNA\tGO:0005634\tnucleus\tcellular_component\n[success]\n")
    monkeypatch.setattr(save_biomart_tables.urllib.request, "urlopen", fake_urlopen)
    output = tmp_path / "biomart_table.tsv"
    monkeypatch.setattr(sys, "argv", ["save_biomart_tables.py", "--biomart_host", "www.ensembl.org",
                                      "--biomart_dataset", "hsapiens_gene_ensembl", "--output", str(output), *options])
    save_biomart_tables.main()
    return output


def test_main_warns_when_biomart_is_another_release_than_the_gtf(monkeypatch, tmp_path, caplog):
    """www.ensembl.org serves the newest release: with the files of an older one the table is written, with a warning."""
    with caplog.at_level("WARNING"):
        output = run_main(monkeypatch, tmp_path, 116, "--ensembl_gtf_url", "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.chr.gtf.gz")
    assert len(output.read_text().splitlines()) == 2
    assert [r.message for r in caplog.records if r.levelname == "WARNING"] == [save_biomart_tables.release_mismatch(116, 115)]


@pytest.mark.parametrize("options", [
    ["--ensembl_gtf_url", "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.chr.gtf.gz"],
    ["--ensembl_gtf_url", "/home/user/my_annotation.gtf"],
    [],
])
def test_main_does_not_warn_when_the_releases_agree_or_are_unknown(monkeypatch, tmp_path, caplog, options):
    with caplog.at_level("WARNING"):
        output = run_main(monkeypatch, tmp_path, 116, *options)
    assert output.exists()
    assert [r.message for r in caplog.records if r.levelname == "WARNING"] == []


def test_main_says_how_to_go_on_when_biomart_fails(monkeypatch, tmp_path, caplog):
    """A BioMart that keeps failing stops the run with the two ways to go on, and leaves no table behind."""
    def fake_urlopen(url, timeout):
        raise OSError("Service unavailable")
    monkeypatch.setattr(save_biomart_tables.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr("time.sleep", lambda seconds: None)
    output = tmp_path / "biomart_table.tsv"
    monkeypatch.setattr(sys, "argv", ["save_biomart_tables.py", "--biomart_host", "e115.ensembl.org",
                                      "--biomart_dataset", "hsapiens_gene_ensembl", "--output", str(output)])
    with pytest.raises(SystemExit) as exit_info:
        save_biomart_tables.main()
    assert exit_info.value.code != 0
    assert "biomart_table=/abs/path/table.tsv" in caplog.text and "layer_go=no reactome=no" in caplog.text
    assert list(tmp_path.iterdir()) == []
