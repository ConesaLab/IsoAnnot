import csv
import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

CODING = ["ENST00000000001", "ENST00000000002", "ENST00000000003", "ENST00000000004", "ENST00000000006"]


def read_tsv(path):
    with open(path) as handle:
        return [line.rstrip("\n").split("\t") for line in handle]


def make_layers(run_script, synthetic_dir, classification, outdir):
    """Exon, junction and NMD layers made by the layer scripts on the fixture, plus a small GO layer."""
    layers = {name: outdir / f"layer_{name}.gtf" for name in ("exons", "junctions", "nmd", "go")}
    for args in (("layer_exons.py", "--gtf_file", synthetic_dir / "ensembl.gtf"),
                 ("layer_junctions.py", "--junctions_file", synthetic_dir / "ensembl_sqanti_junctions.txt",
                  "--classification_file", classification),
                 ("layer_nmd.py", "--nmd_file", synthetic_dir / "ensembl_sqanti_NMD.txt",
                  "--classification_file", classification)):
        name = args[0][len("layer_"):-len(".py")]
        result = run_script(*args, "--output", layers[name])
        assert result.returncode == 0, result.stderr
    layers["go"].write_text("ENST00000000001\tGeneOntology\tC\t.\t.\t.\t.\t.\tID=GO:0005737; Name=cytoplasm\n"
                            "ENST00000000003\tGeneOntology\tF\t.\t.\t.\t.\t.\tID=GO:0005515; Name=protein binding\n")
    return layers


def run_t2go(run_script, synthetic_dir, tmp_path, classification=None, protein_assoc=None):
    """Runs t2goAnnotationFile.py as rule tappas_annotation does and returns the output rows."""
    classification = classification or synthetic_dir / "ensembl_sqanti_classification.txt"
    layers = make_layers(run_script, synthetic_dir, classification, tmp_path)
    output = tmp_path / "annotation.gff3"
    args = ["--classification_file", classification, "--input_transcripts", layers["nmd"],
            "--input_genomic", layers["exons"], layers["junctions"], "--input_protein", layers["go"],
            "--output", output]
    if protein_assoc:
        args += ["--protein_association", protein_assoc]
    result = run_script("t2goAnnotationFile.py", *args)
    assert result.returncode == 0, result.stderr
    return read_tsv(output)


def test_reference_sqanti_orf_length_is_an_integer(synthetic_dir):
    """referenceSQANTI.py writes ORF_length as an integer number of codons, not 150.0 (GitHub #10)."""
    with open(synthetic_dir / "ensembl_sqanti_classification.txt") as handle:
        orf_length = {row["isoform"]: row["ORF_length"] for row in csv.DictReader(handle, delimiter="\t")}
    assert orf_length == {"ENST00000000001": "150", "ENST00000000002": "87", "ENST00000000003": "187",
                          "ENST00000000004": "187", "ENST00000000005": "NA", "ENST00000000006": "267"}


@pytest.mark.parametrize("float_orf_length", [False, True])
def test_protein_end_is_an_integer(synthetic_dir, run_script, tmp_path, float_orf_length):
    """
    The protein end is an integer: IsoAnnotLite crashes on int('180.0') and tappAS reads it as 0 (GitHub
    #10). Classifications made before the fix still have float ORF lengths.
    """
    classification = synthetic_dir / "ensembl_sqanti_classification.txt"
    if float_orf_length:
        with open(classification) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            if row["ORF_length"] != "NA":
                row["ORF_length"] = f"{int(row['ORF_length'])}.0"
        classification = tmp_path / "classification.txt"
        with open(classification, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    rows = run_t2go(run_script, synthetic_dir, tmp_path, classification=classification)
    proteins = {row[0]: (row[3], row[4]) for row in rows if row[2] == "protein"}
    assert proteins == {"ENST00000000001": ("1", "150"), "ENST00000000002": ("1", "87"),
                        "ENST00000000003": ("1", "187"), "ENST00000000004": ("1", "187"),
                        "ENST00000000006": ("1", "267")}
