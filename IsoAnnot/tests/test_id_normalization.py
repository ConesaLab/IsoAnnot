import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from IsoAnnot import normalize_id


@pytest.mark.parametrize("identifier, expected", [
    ("ENST00000000233.5", "ENST00000000233"),
    ("ENSP00000495982.2", "ENSP00000495982"),
    ("ENSMUST00000193812.2", "ENSMUST00000193812"),
    ("FBtr0070000.1", "FBtr0070000"),
    ("NM_000005.3", "NM_000005"),
    ("NP_000005.3", "NP_000005"),
    ("XM_011541469.2", "XM_011541469"),
    ("NR_046018.2", "NR_046018"),
    ("ENST00000000233", "ENST00000000233"),
    ("NM_000005", "NM_000005"),
    ("PB.1.1", "PB.1.1"),
    ("PB.12.3", "PB.12.3"),
    ("STRG.9.1", "STRG.9.1"),
    ("AT1G01010.1", "AT1G01010.1"),
    ("transcript_1", "transcript_1"),
])
def test_normalize_id(identifier, expected):
    """Only Ensembl, FlyBase and RefSeq versions are stripped; PacBio-style IDs keep their dots."""
    assert normalize_id(identifier) == expected


def write_nls_file(path, ids):
    """NucImport predictions as expand_nls leaves them (nls_final_expanded.tsv)."""
    lines = ["ID\tImpProb\tcNLSProb\tSeq\tClass\tPos"]
    lines += [f"{prot_id}\t0.900\t0.800\tKRKRKRKRKRKRKRKRKRKR\tmonopartite\t12" for prot_id in ids]
    path.write_text("\n".join(lines) + "\n")


def nls_row(isoform):
    return f"{isoform}\tNucImport\tNLS\t12\t32\t.\t.\t.\tID=NLS_monopartite; Name=NLS_monopartite; Desc=NLS_NucImport"


def test_layer_nls_versioned_ensembl_proteins(tmp_path, synthetic_dir, run_script):
    """pep.all IDs are versioned and refProt is not: the NLS layer was empty in ensembl mode."""
    nls_file = tmp_path / "nls_final_expanded.tsv"
    write_nls_file(nls_file, ["ENSP00000000001.2", "ENSP00000000006.3"])
    output = tmp_path / "layer_nls.gtf"
    result = run_script("layer_nls.py", "--nls_file", nls_file,
                        "--classification_file", synthetic_dir / "ensembl_sqanti_classification.txt",
                        "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [nls_row("ENST00000000001"), nls_row("ENST00000000006")]


def test_layer_nls_keeps_pacbio_ids(tmp_path, run_script):
    """The version fallback must not cut PacBio ORF IDs (PB.1.1) at the first dot."""
    classification = tmp_path / "classification.txt"
    classification.write_text("isoform\tcoding\nPB.1.1\tcoding\nPB.12.3\tcoding\n")
    nls_file = tmp_path / "nls_final_expanded.tsv"
    write_nls_file(nls_file, ["PB.1.1", "PB.12.3"])
    output = tmp_path / "layer_nls.gtf"
    result = run_script("layer_nls.py", "--nls_file", nls_file, "--classification_file", classification,
                        "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [nls_row("PB.1.1"), nls_row("PB.12.3")]
