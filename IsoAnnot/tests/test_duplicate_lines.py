import os
import sys
import types
from unittest.mock import MagicMock

import pandas as pd
import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)


def test_layer_interproscan_writes_each_clan_once(run_script, synthetic_dir, tmp_path):
    """Several domains of one clan give one domain line each but one clan line per transcript."""
    interproscan = tmp_path / "interproResultsPfam.tsv"
    interproscan.write_text(
        "ENSP00000000001\tPfam\t-\tPF13855\tLRR_8\tLeucine rich repeat\t10\t40\n"
        "ENSP00000000001\tPfam\t-\tPF13855\tLRR_8\tLeucine rich repeat\t60\t90\n"
        "ENSP00000000002\tPfam\t-\tPF13855\tLRR_8\tLeucine rich repeat\t5\t30\n"
        "ENSP00000000002\tPfam\t-\tPF12799\tLRR_4\tLeucine Rich repeats (2 copies)\t40\t70\n"
    )
    pfam = tmp_path / "Pfam-A.clans.tsv"
    pfam.write_text("PF13855\tCL0022\tLRR\tLRR_8\tLeucine rich repeat\n"
                    "PF12799\tCL0022\tLRR\tLRR_4\tLeucine Rich repeats (2 copies)\n")
    clans = tmp_path / "Pfam-C"
    clans.write_text("# STOCKHOLM 1.0\n#=GF ID   LRR\n#=GF AC   CL0022.41\n#=GF DE   Leucine Rich Repeat\n//\n")
    output = tmp_path / "layer_interproscan.gtf"
    result = run_script("layer_interproscan.py", "--interproscan_file", interproscan, "--pfam_file", pfam, "--pfam_clan_desc", clans,
                        "--classification_file", synthetic_dir / "ensembl_sqanti_classification.txt",
                        "--output", output)
    assert result.returncode == 0, result.stderr
    domain = "PFAM\tDOMAIN\t{}\t{}\t.\t.\t.\tID=PF13855; Name=LRR_8; Desc=Leucine rich repeat"
    clan = "PFAM\tclan\t.\t.\t.\t.\t.\tID=CL0022; Name=LRR; Desc=Leucine Rich Repeat"
    # LRR_4 is a second family of clan CL0022: ENST00000000002 already has that clan line
    assert output.read_text().splitlines() == [
        "ENST00000000001\t" + domain.format(10, 40),
        "ENST00000000001\t" + clan,
        "ENST00000000001\t" + domain.format(60, 90),
        "ENST00000000002\t" + domain.format(5, 30),
        "ENST00000000002\t" + clan,
        "ENST00000000002\tPFAM\tDOMAIN\t40\t70\t.\t.\t.\tID=PF12799; Name=LRR_4; Desc=Leucine Rich repeats (2 copies)",
    ]


def run_layer_reactome(monkeypatch, tmp_path, classification, reactome_text, gene_names, db):
    """
    Runs layer_reactome.main() in this process, with BioMart's gene ID-name answer given by the test:
    pybiomart is not needed and nothing goes to the network. Returns the output lines.
    """
    pybiomart = types.ModuleType("pybiomart")
    pybiomart.Dataset = MagicMock()
    monkeypatch.setitem(sys.modules, "pybiomart", pybiomart)
    sys.modules.pop("layer_reactome", None)
    import layer_reactome
    monkeypatch.setattr(layer_reactome, "query_biomart_with_retry", lambda **kwargs: gene_names.copy())
    reactome = tmp_path / "Ensembl2Reactome_All_Levels.txt"
    reactome.write_text(reactome_text)
    output = tmp_path / "layer_reactome.gtf"
    monkeypatch.setattr(sys, "argv", ["layer_reactome.py", "--reactome_file", str(reactome), "--classification_file", str(classification),
                                      "--biomart_host", "www.ensembl.org", "--biomart_dataset", "hsapiens_gene_ensembl",
                                      "--species", "Homo sapiens", "--db", db, "--output", str(output)])
    with pytest.raises(SystemExit) as exit_info:
        layer_reactome.main()
    assert exit_info.value.code == 0
    return output.read_text().splitlines()


# R-HSA-1 is listed twice for ENSG00000000001 (TAS and IEA evidence) and once more for ENSG00000000099,
# a second Ensembl ID with the gene name SYNA
REACTOME = (
    "ENSG00000000001\tR-HSA-1\thttps://reactome.org/R-HSA-1\tPathway one\tTAS\tHomo sapiens\n"
    "ENSG00000000001\tR-HSA-1\thttps://reactome.org/R-HSA-1\tPathway one\tIEA\tHomo sapiens\n"
    "ENSG00000000099\tR-HSA-1\thttps://reactome.org/R-HSA-1\tPathway one\tTAS\tHomo sapiens\n"
    "ENSG00000000001\tR-HSA-2\thttps://reactome.org/R-HSA-2\tPathway two\tTAS\tHomo sapiens\n"
)
GENE_NAMES = pd.DataFrame({"Gene stable ID": ["ENSG00000000001", "ENSG00000000099"], "Gene name": ["SYNA", "SYNA"]})


def test_layer_reactome_writes_each_pathway_once(monkeypatch, synthetic_dir, tmp_path):
    """A pathway listed once per evidence code, or under two Ensembl IDs of one gene name, is written once per transcript."""
    lines = run_layer_reactome(monkeypatch, tmp_path, synthetic_dir / "ensembl_sqanti_classification.txt", REACTOME, GENE_NAMES, "ensembl")
    assert lines == [
        f"{transcript}\tREACTOME\tPATHWAY\t.\t.\t.\t.\t.\tID={pathway}; Name={name}"
        for transcript in ["ENST00000000001", "ENST00000000002"]
        for pathway, name in [("R-HSA-1", "Pathway one"), ("R-HSA-2", "Pathway two")]
    ]


def test_layer_reactome_mytranscripts_each_pathway_once(monkeypatch, synthetic_dir, tmp_path):
    """In mytranscripts mode genes are matched by Ensembl ID; the evidence-code repeats are dropped there too."""
    classification = tmp_path / "classification.txt"
    classification.write_text((synthetic_dir / "ensembl_sqanti_classification.txt").read_text().replace("\tSYNA\t", "\tENSG00000000001\t"))
    lines = run_layer_reactome(monkeypatch, tmp_path, classification, REACTOME, GENE_NAMES, "mytranscripts")
    assert lines == [
        f"{transcript}\tREACTOME\tPATHWAY\t.\t.\t.\t.\t.\tID={pathway}; Name={name}"
        for transcript in ["ENST00000000001", "ENST00000000002"]
        for pathway, name in [("R-HSA-1", "Pathway one"), ("R-HSA-2", "Pathway two")]
    ]
