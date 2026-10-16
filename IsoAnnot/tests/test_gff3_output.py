import csv
import json
import os
import re
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from IsoAnnot import natural_key
from t2goAnnotationFile import resolve_protein_id
from transcript2reference import get_novel_protein_ids

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


def write_classification(path, rows):
    """Writes a classification with the columns t2goAnnotationFile.py and reference_protein_assoc.py read."""
    columns = ["isoform", "chrom", "strand", "length", "structural_category", "associated_gene",
               "associated_transcript", "coding", "ORF_length", "CDS_start", "CDS_end", "refProt"]
    with open(path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(v) for v in row) + "\n")
    return path


def run_reference_protein_assoc(run_script, synthetic_dir, tmp_path, database, classification=None):
    classification = classification or synthetic_dir / "ensembl_sqanti_classification.txt"
    output = tmp_path / "reference_protein_assoc_data.txt"
    result = run_script("reference_protein_assoc.py", "--classification_file", classification,
                        "--species_db", synthetic_dir / "parsed_databases.json", "--database", database,
                        "--output", output)
    assert result.returncode == 0, result.stderr
    return output, result


def test_natural_key():
    """Numbers inside IDs compare as numbers (GitHub #10)."""
    ids = ["Q9HBE1-10", "NP_000005.10", "Q9HBE1-2", "NP_000005.3", "Q9HBE1", "PB.10.1", "PB.2.1"]
    assert sorted(ids, key=natural_key) == ["NP_000005.3", "NP_000005.10", "PB.2.1", "PB.10.1",
                                            "Q9HBE1", "Q9HBE1-2", "Q9HBE1-10"]


@pytest.mark.parametrize("proteins,expected", [
    # SwissProt beats TrEMBL, which beats the native IDs
    ([("ENSP00000000001", "ensembl"), ("A0A024R1K6", "trembl"), ("Q9BZE7", "swissprot")],
     ("Q9BZE7", ["A0A024R1K6", "ENSP00000000001"])),
    # A canonical entry beats an isoform; isoform numbers compare as numbers
    ([("Q9HBE1-10", "swissprot"), ("Q9HBE1-2", "swissprot"), ("Q9HBE1", "swissprot")],
     ("Q9HBE1", ["Q9HBE1-2", "Q9HBE1-10"])),
    ([("Q9HBE1-10", "swissprot"), ("Q9HBE1-2", "swissprot")], ("Q9HBE1-2", ["Q9HBE1-10"])),
    # Any canonical entry beats an isoform entry, even one with a lower accession
    ([("P12344-2", "swissprot"), ("P12345", "swissprot")], ("P12345", ["P12344-2"])),
    # Without UniProt, the lowest native ID, whatever the transcript's own protein is
    ([("ENSP00000384194", "ensembl"), ("ENSP00000350199", "ensembl")], ("ENSP00000350199", ["ENSP00000384194"])),
    # mytranscripts: Ensembl before RefSeq, and Protein_<isoform> only when nothing else matched
    ([("NP_000005", "refseq"), ("ENSP00000325785", "ensembl")], ("ENSP00000325785", ["NP_000005"])),
    ([("Protein_PB.1.1", "novel")], ("Protein_PB.1.1", [])),
    # Duplicated rows are reported once
    ([("P0SYA1", "swissprot"), ("P0SYA1", "swissprot")], ("P0SYA1", [])),
])
def test_resolve_protein_id(proteins, expected):
    """The protein ID and Desc only depend on the set of same-sequence IDs (GitHub #10)."""
    assert resolve_protein_id(proteins) == expected
    assert resolve_protein_id(list(reversed(proteins))) == expected


def test_reference_protein_assoc_ensembl(synthetic_dir, run_script, tmp_path):
    """Ensembl mode links each transcript to the UniProt and Ensembl IDs of its protein sequence."""
    output, _ = run_reference_protein_assoc(run_script, synthetic_dir, tmp_path, "ensembl")
    assert read_tsv(output) == [
        ["ENST00000000001", "P0SYA1", "swissprot", "sequence"],
        ["ENST00000000001", "ENSP00000000001", "ensembl", "sequence"],
        ["ENST00000000002", "A0A0S1YNA6", "trembl", "sequence"],
        ["ENST00000000002", "ENSP00000000002", "ensembl", "sequence"],
        ["ENST00000000003", "P0SYB1-2", "swissprot", "sequence"],
        ["ENST00000000003", "ENSP00000000003", "ensembl", "sequence"],
        ["ENST00000000003", "ENSP00000000004", "ensembl", "sequence"],
        ["ENST00000000004", "P0SYB1-2", "swissprot", "sequence"],
        ["ENST00000000004", "ENSP00000000003", "ensembl", "sequence"],
        ["ENST00000000004", "ENSP00000000004", "ensembl", "sequence"],
        ["ENST00000000006", "ENSP00000000006", "ensembl", "sequence"],
    ]


def test_reference_protein_assoc_refseq(synthetic_dir, run_script, tmp_path):
    """RefSeq mode looks up versioned refProt IDs and keeps RefSeq, not Ensembl, IDs."""
    classification = write_classification(tmp_path / "classification.txt", [
        ["NM_900003.2", "1", "-", 710, "protein_coding", "SYNB", "NM_900003.2", "coding", 187, 51, 611, "NP_900003.1"],
        ["NR_900005.1", "1", "+", 700, "lncRNA", "SYNC", "NR_900005.1", "non_coding", "NA", "NA", "NA", "NA"],
        ["NM_900006.3", "1", "+", 1100, "protein_coding", "SYND", "NM_900006.3", "coding", 267, 101, 901, "NP_900006.2"],
        ["NM_900009.1", "1", "+", 900, "protein_coding", "SYNE", "NM_900009.1", "coding", 200, 101, 700, "NP_900009.1"],
    ])
    output, result = run_reference_protein_assoc(run_script, synthetic_dir, tmp_path, "refseq", classification)
    assert read_tsv(output) == [
        ["NM_900003.2", "P0SYB1-2", "swissprot", "sequence"],
        ["NM_900003.2", "NP_900003", "refseq", "sequence"],
        ["NM_900003.2", "NP_900004", "refseq", "sequence"],
        ["NM_900006.3", "NP_900006", "refseq", "sequence"],
    ]
    assert "1 of 3 coding transcripts" in result.stderr


def test_reference_protein_assoc_ambiguous_native_id(run_script, tmp_path):
    """
    A native ID that maps to two sequences once the version is removed (AT1G01010.1 and AT1G01010.2 in
    Ensembl Plants) is not used to pick a sequence.
    """
    species_db = tmp_path / "parsed_databases.json"
    species_db.write_text(json.dumps({"MAAA": {"swissprot": ["Q00001"], "ensembl": ["AT1G01010"]},
                                      "MCCC": {"ensembl": ["AT1G01010"]},
                                      "MGGG": {"swissprot": ["Q00002"], "ensembl": ["AT1G01020"]}}))
    classification = write_classification(tmp_path / "classification.txt", [
        ["AT1G01010.1", "1", "+", 900, "protein_coding", "G1", "AT1G01010.1", "coding", 3, 1, 9, "AT1G01010.1"],
        ["AT1G01020.1", "1", "+", 900, "protein_coding", "G2", "AT1G01020.1", "coding", 3, 1, 9, "AT1G01020.1"],
    ])
    output = tmp_path / "assoc.txt"
    result = run_script("reference_protein_assoc.py", "--classification_file", classification,
                        "--species_db", species_db, "--database", "ensembl", "--output", output)
    assert result.returncode == 0, result.stderr
    assert read_tsv(output) == [["AT1G01020.1", "Q00002", "swissprot", "sequence"],
                                ["AT1G01020.1", "AT1G01020", "ensembl", "sequence"]]


def test_novel_protein_ids():
    """ORFs with no exact match get Protein_<isoform>, shared by isoforms with the same ORF."""
    isoform_proteins = {"PB.10.1": "MKKK", "PB.2.1": "MKKK", "PB.3.1": "MAAA", "PB.4.1": "MLLL",
                        "PB:5,1": "MWWW"}
    assert get_novel_protein_ids(isoform_proteins, {"MAAA": {"swissprot": ["Q00001"]}}) == {
        "PB.10.1": "Protein_PB.2.1", "PB.2.1": "Protein_PB.2.1", "PB.4.1": "Protein_PB.4.1",
        "PB:5,1": "Protein_PB_5_1"}


def test_transcript2reference_writes_novel_rows(synthetic_dir, run_script, tmp_path):
    """transcript2reference.py writes the shared Protein_<isoform> ID of novel ORFs to protein_assoc."""
    with open(synthetic_dir / "parsed_databases.json") as handle:
        p0sya1 = next(seq for seq, dbs in json.load(handle).items() if dbs.get("swissprot") == ["P0SYA1"])
    corrected_gtf = tmp_path / "corrected.gtf"
    corrected_gtf.write_text("".join(f'1\tPacBio\texon\t1001\t2000\t.\t+\t.\ttranscript_id "{isoform}"; gene_id "PB.1";\n'
                                     for isoform in ("PB.2.1", "PB.3.1", "PB.10.1")))
    orfs = tmp_path / "corrected.faa"
    orfs.write_text(f">PB.2.1 GeneMark.hmm|4_aa|+|11|22\nMKKK\n"
                    f">PB.3.1 GeneMark.hmm|150_aa|+|11|460\n{p0sya1}\n"
                    f">PB.10.1 GeneMark.hmm|4_aa|+|11|22\nMKKK\n")
    protein_assoc = tmp_path / "protein_assoc_data.txt"
    nmd = tmp_path / "nmd.txt"
    result = run_script("transcript2reference.py", "--ensembl_gtf", synthetic_dir / "ensembl.gtf",
                        "--corrected_gtf", corrected_gtf, "--corrected_fasta_proteins", orfs,
                        "--database", "mytranscripts", "--species_db", synthetic_dir / "parsed_databases.json",
                        "--output_assoc", protein_assoc, "--output_nmd", nmd)
    assert result.returncode == 0, result.stderr
    assert [line.split("\t") for line in protein_assoc.read_text().splitlines()] == [
        ["PB.2.1", "Protein_PB.2.1", "novel", "sequence"],
        ["PB.3.1", "P0SYA1", "swissprot", "sequence"],
        ["PB.3.1", "NP_900001", "refseq", "sequence"],
        ["PB.3.1", "ENSP00000000001", "ensembl", "sequence"],
        ["PB.10.1", "Protein_PB.2.1", "novel", "sequence"],
    ]
    assert sorted(os.listdir(tmp_path)) == ["corrected.faa", "corrected.gtf", "nmd.txt", "protein_assoc_data.txt"]


def test_transcript2reference_orfs_with_stop_codon(synthetic_dir, run_script, tmp_path):
    """SQANTI3 ORFs end with *, the reference proteins don't: they must still match, so long-read
    proteins get their UniProt/Ensembl/RefSeq IDs, the same as for an ORF written without the *."""
    with open(synthetic_dir / "parsed_databases.json") as handle:
        p0sya1 = next(seq for seq, dbs in json.load(handle).items() if dbs.get("swissprot") == ["P0SYA1"])
    corrected_gtf = tmp_path / "corrected.gtf"
    corrected_gtf.write_text('1\tPacBio\texon\t1001\t2000\t.\t+\t.\ttranscript_id "PB.3.1"; gene_id "PB.1";\n')
    outputs = {}
    for name, stop in (("with_stop", "*"), ("without_stop", "")):
        orfs = tmp_path / f"{name}.faa"
        orfs.write_text(f">PB.3.1 GeneMark.hmm|150_aa|+|11|460\n{p0sya1}{stop}\n")
        result = run_script("transcript2reference.py", "--ensembl_gtf", synthetic_dir / "ensembl.gtf",
                            "--corrected_gtf", corrected_gtf, "--corrected_fasta_proteins", orfs,
                            "--database", "mytranscripts", "--species_db", synthetic_dir / "parsed_databases.json",
                            "--output_assoc", tmp_path / f"{name}_assoc.txt", "--output_nmd", tmp_path / f"{name}_nmd.txt")
        assert result.returncode == 0, result.stderr
        outputs[name] = (tmp_path / f"{name}_assoc.txt").read_text()
    assert [line.split("\t") for line in outputs["with_stop"].splitlines()] == [
        ["PB.3.1", "P0SYA1", "swissprot", "sequence"],
        ["PB.3.1", "NP_900001", "refseq", "sequence"],
        ["PB.3.1", "ENSP00000000001", "ensembl", "sequence"],
    ]
    assert outputs["with_stop"] == outputs["without_stop"]


def test_full_protein_ids_and_desc(synthetic_dir, run_script, tmp_path):
    """
    CDS and protein lines carry the full UniProt ID, with the other same-sequence IDs in Desc, instead of
    the first character of refProt (GitHub #10). SYNB-201 and SYNB-202 share one ID and Desc.
    """
    protein_assoc, _ = run_reference_protein_assoc(run_script, synthetic_dir, tmp_path, "ensembl")
    rows = run_t2go(run_script, synthetic_dir, tmp_path, protein_assoc=protein_assoc)
    cds = {row[0]: row[8] for row in rows if row[2] == "CDS"}
    proteins = {row[0]: row[8] for row in rows if row[2] == "protein"}
    assert cds == {
        "ENST00000000001": "ID=P0SYA1; Desc=ENSP00000000001",
        "ENST00000000002": "ID=A0A0S1YNA6; Desc=ENSP00000000002",
        "ENST00000000003": "ID=P0SYB1-2; Desc=ENSP00000000003,ENSP00000000004",
        "ENST00000000004": "ID=P0SYB1-2; Desc=ENSP00000000003,ENSP00000000004",
        "ENST00000000006": "ID=ENSP00000000006",
    }
    assert proteins == cds


TAPPAS_FEATURES = {"transcript": "t", "gene": "g", "CDS": "c", "genomic": "G", "exon": "e",
                   "splice_junction": "j", "protein": "p"}


def block_letter(row):
    """One letter per line: tappAS features by name, layer features by their PosType (T, P or N)."""
    return TAPPAS_FEATURES[row[2]] if row[1] == "tappAS" else row[8][-1]


def test_end_to_end_gff3(synthetic_dir, run_script, tmp_path):
    """
    t2goAnnotationFile.py and renameFeatures.py on the fixture, as rules tappas_annotation and renameFeatures
    run them, give the layout of the legacy tappAS files: each transcript starts with its transcript line
    (tappAS needs it first),
    then the transcript, genomic and protein blocks; full protein IDs and integer coordinates
    (GitHub #10).
    """
    protein_assoc, _ = run_reference_protein_assoc(run_script, synthetic_dir, tmp_path, "ensembl")
    run_t2go(run_script, synthetic_dir, tmp_path, protein_assoc=protein_assoc)
    result = run_script("renameFeatures.py", tmp_path / "annotation.gff3")
    assert result.returncode == 0, result.stderr
    rows = read_tsv(tmp_path / "annotation.gff3_mod")

    assert all(len(row) == 9 and re.search(r"[^;\s]; PosType=[TGPN]$", row[8]) for row in rows)

    # Each transcript's lines are contiguous and in the legacy tappAS order
    transcripts = list(dict.fromkeys(row[0] for row in rows))
    assert transcripts == CODING[:4] + ["ENST00000000005", "ENST00000000006"]
    assert [row[0] for row in rows] == sorted((row[0] for row in rows), key=transcripts.index)
    blocks = {transcript: "".join(block_letter(row) for row in rows if row[0] == transcript)
              for transcript in transcripts}
    assert blocks == {"ENST00000000001": "tgcGeeejjpN", "ENST00000000002": "tgcTGeeeejjjp",
                      "ENST00000000003": "tgcGeeeejjjpN", "ENST00000000004": "tgcGeeeeejjjjp",
                      "ENST00000000005": "tgGeej", "ENST00000000006": "tgcGep"}
    assert all(re.match(r"^tgc?T*Ge+j*(p[NP]*)?$", block) for block in blocks.values())

    # Integer coordinates, except the position-less GO and NMD features
    for row in rows:
        if row[1] in ("GeneOntology", "NMD"):
            assert (row[3], row[4]) == (".", ".")
        else:
            assert re.match(r"^\d+$", row[3]) and re.match(r"^\d+$", row[4]), row

    # Same ID and Desc on CDS and protein lines, shared by the two SYNB transcripts
    cds = {row[0]: (row[3], row[4], row[8]) for row in rows if row[2] == "CDS"}
    proteins = {row[0]: (row[3], row[4], row[8]) for row in rows if row[2] == "protein"}
    assert cds == {
        "ENST00000000001": ("101", "550", "ID=P0SYA1; Desc=ENSP00000000001; PosType=T"),
        "ENST00000000002": ("101", "361", "ID=A0A0S1YNA6; Desc=ENSP00000000002; PosType=T"),
        "ENST00000000003": ("51", "611", "ID=P0SYB1-2; Desc=ENSP00000000003,ENSP00000000004; PosType=T"),
        "ENST00000000004": ("111", "671", "ID=P0SYB1-2; Desc=ENSP00000000003,ENSP00000000004; PosType=T"),
        "ENST00000000006": ("101", "901", "ID=ENSP00000000006; PosType=T"),
    }
    assert proteins == {
        "ENST00000000001": ("1", "150", "ID=P0SYA1; Desc=ENSP00000000001; PosType=P"),
        "ENST00000000002": ("1", "87", "ID=A0A0S1YNA6; Desc=ENSP00000000002; PosType=P"),
        "ENST00000000003": ("1", "187", "ID=P0SYB1-2; Desc=ENSP00000000003,ENSP00000000004; PosType=P"),
        "ENST00000000004": ("1", "187", "ID=P0SYB1-2; Desc=ENSP00000000003,ENSP00000000004; PosType=P"),
        "ENST00000000006": ("1", "267", "ID=ENSP00000000006; PosType=P"),
    }
    assert all(re.match(r"^ID=[A-Za-z0-9_.\-]+;", attributes) for _, _, attributes in proteins.values())


def test_t2go_empty_blocks(run_script, tmp_path):
    """
    Layers can be empty (no UTR motif, no NMD transcript). When every file of a block is empty, t2go
    writes the primary lines instead of crashing on pd.concat([]).
    """
    classification = write_classification(tmp_path / "classification.txt", [
        ["ENST00000000005", "1", "+", 700, "lncRNA", "SYNC", "ENST00000000005", "non_coding", "NA", "NA", "NA", "NA"],
    ])
    empty = tmp_path / "empty.gtf"
    empty.write_text("")
    exons = tmp_path / "exons.gtf"
    exons.write_text("ENST00000000005\ttappAS\texon\t12001\t12400\t.\t+\t.\tChr=chr1\n")
    output = tmp_path / "annotation.gff3"
    result = run_script("t2goAnnotationFile.py", "--classification_file", classification,
                        "--input_transcripts", empty, empty, "--input_genomic", exons, "--input_protein", empty,
                        "--output", output)
    assert result.returncode == 0, result.stderr
    assert [row[:5] for row in read_tsv(output)] == [
        ["ENST00000000005", "tappAS", "transcript", "1", "700"],
        ["ENST00000000005", "tappAS", "gene", "1", "700"],
        ["ENST00000000005", "tappAS", "genomic", "1", "1"],
        ["ENST00000000005", "tappAS", "exon", "12001", "12400"],
    ]


def test_t2go_warns_about_dropped_layer_lines(run_script, tmp_path):
    """
    Layer lines whose transcript ID is not in the classification (e.g. versioned vs unversioned IDs) are
    dropped. t2go now says how many, once per layer, instead of emptying the layer silently.
    """
    classification = write_classification(tmp_path / "classification.txt", [
        ["ENST00000000005", "1", "+", 700, "lncRNA", "SYNC", "ENST00000000005", "non_coding", "NA", "NA", "NA", "NA"],
    ])
    exons = tmp_path / "layer_exons.gtf"
    exons.write_text("ENST00000000005\ttappAS\texon\t12001\t12400\t.\t+\t.\tChr=chr1\n"
                     "ENST00000000005\ttappAS\texon\t13001\t13300\t.\t+\t.\tChr=chr1\n")
    repeats = tmp_path / "layer_repeatmasker.gtf"
    repeats.write_text("ENST00000000005.2\tRepeatMasker\trepeat\t10\t60\t.\t+\t.\tID=AluY\n"
                       "ENST00000000005.2\tRepeatMasker\trepeat\t90\t120\t.\t+\t.\tID=L1\n"
                       "ENST00000000005\tRepeatMasker\trepeat\t200\t260\t.\t+\t.\tID=MIR\n")
    output = tmp_path / "annotation.gff3"
    result = run_script("t2goAnnotationFile.py", "--classification_file", classification,
                        "--input_transcripts", repeats, "--input_genomic", exons, "--input_protein", exons,
                        "--output", output)
    assert result.returncode == 0, result.stderr
    assert [row[8] for row in read_tsv(output) if row[1] == "RepeatMasker"] == ["ID=MIR"]
    warnings = [line for line in result.stderr.splitlines() if "not in the classification file" in line]
    assert warnings == [f"WARNING:root:2 of 3 lines in {repeats} were dropped because their transcript ID is "
                        f"not in the classification file. If that is most of the file, check the "
                        f"transcript_versioned setting of the config file."]


def test_rename_features_keeps_unknown_pairs(run_script, tmp_path):
    """
    A source/feature pair renameFeatures.py does not know (e.g. from a user protein_gtf) no longer stops
    the loop and truncates the file with exit 0. The line is kept as PosType=N with one summary warning.
    """
    gff3 = tmp_path / "annotation.gff3"
    gff3.write_text("TX1\ttappAS\ttranscript\t1\t750\t.\t+\t.\tID=TX1\n"
                    "TX1\tMyTool\tmotif\t10\t20\t.\t+\t.\tID=M1\n"
                    "TX1\tPFAM\tFAMILY\t1\t50\t.\t.\t.\tID=PF00001\n"
                    "TX1\tMyTool\tmotif\t30\t40\t.\t+\t.\tID=M2\n"
                    "TX1\ttappAS\tprotein\t1\t150\t.\t+\t.\tID=P0SYA1\n")
    result = run_script("renameFeatures.py", gff3)
    assert result.returncode == 0, result.stderr
    assert [row[8] for row in read_tsv(tmp_path / "annotation.gff3_mod")] == [
        "ID=TX1; PosType=T", "ID=M1; PosType=N", "ID=PF00001; PosType=N", "ID=M2; PosType=N",
        "ID=P0SYA1; PosType=P"]
    warnings = [line for line in result.stderr.splitlines() if "renameFeatures.py does not know" in line]
    assert warnings == ["WARNING:root:3 lines have a source/feature pair that renameFeatures.py does not know "
                        "(MyTool/motif: 2, PFAM/FAMILY: 1). They were kept with PosType=N (no position)."]


@pytest.mark.parametrize("feature", ["CA_BIND", "NONT_STD", "NON_STD", "ACT_SITE"])
def test_rename_features_uniprot_features(run_script, tmp_path, feature):
    """Every feature layer_uniprot can write is a positional protein feature."""
    gff3 = tmp_path / "annotation.gff3"
    gff3.write_text(f"TX1\tUniProtKB/Swiss-Prot_Phosphosite\t{feature}\t5\t9\t.\t.\t.\tID=X\n")
    result = run_script("renameFeatures.py", gff3)
    assert result.returncode == 0, result.stderr
    assert read_tsv(tmp_path / "annotation.gff3_mod") == [
        ["TX1", "UniProtKB/Swiss-Prot_Phosphosite", feature, "5", "9", ".", ".", ".", "ID=X; PosType=P"]]
    assert "does not know" not in result.stderr


def test_rename_features_malformed_line(run_script, tmp_path):
    """
    A line without 9 columns stops with its line number instead of a NameError on idx, and leaves no
    partial output that a later run would reuse.
    """
    gff3 = tmp_path / "annotation.gff3"
    gff3.write_text("TX1\ttappAS\ttranscript\t1\t750\t.\t+\t.\tID=TX1\n"
                    "TX1\ttappAS\tgene\t1\t750\n")
    result = run_script("renameFeatures.py", gff3)
    assert result.returncode != 0
    assert "Invalid line 2 (5 columns instead of 9)" in result.stderr
    assert "NameError" not in result.stderr
    assert sorted(os.listdir(tmp_path)) == ["annotation.gff3"]


def test_protein_id_fallbacks(run_script, tmp_path):
    """
    Without a sequence match the ID is the refProt (reference modes) or Protein_<isoform> (mytranscripts),
    never a single character or a shared "Novel" (GitHub #10).
    """
    classification = write_classification(tmp_path / "classification.txt", [
        ["ENST00000000001", "1", "+", 750, "protein_coding", "SYNA", "ENST00000000001", "coding", 150, 101, 550, "ENSP00000000001"],
        ["ENST00000000006", "1", "+", 1100, "protein_coding", "SYND", "ENST00000000006", "coding", 267, 101, 901, "ENSP00000000006"],
    ])
    protein_assoc = tmp_path / "assoc.txt"
    protein_assoc.write_text("ENST00000000001\tP0SYA1\tswissprot\tsequence\n")
    exons = tmp_path / "exons.gtf"
    exons.write_text("ENST00000000001\ttappAS\texon\t1001\t1750\t.\t+\t.\tChr=chr1\n")
    output = tmp_path / "annotation.gff3"
    result = run_script("t2goAnnotationFile.py", "--classification_file", classification,
                        "--input_transcripts", exons, "--input_genomic", exons, "--input_protein", exons,
                        "--protein_association", protein_assoc, "--output", output)
    assert result.returncode == 0, result.stderr
    cds = {row[0]: row[8] for row in read_tsv(output) if row[2] == "CDS"}
    assert cds == {"ENST00000000001": "ID=P0SYA1", "ENST00000000006": "ID=ENSP00000000006"}
    assert ("WARNING:root:1 coding transcripts have no reference protein with the same sequence. "
            "Their refProt from the classification file is used as protein ID.") in result.stderr.splitlines()

    # mytranscripts: SQANTI3 classifications have no refProt column
    sqanti = tmp_path / "sqanti_classification.txt"
    sqanti.write_text("isoform\tchrom\tstrand\tlength\tstructural_category\tassociated_gene\tassociated_transcript"
                      "\tcoding\tORF_length\tCDS_start\tCDS_end\n"
                      "PB.1.1\t1\t+\t750\tnovel_not_in_catalog\tSYNA\tnovel\tcoding\t150\t101\t553\n"
                      "PB:2,1\t1\t+\t750\tnovel_not_in_catalog\tSYNA\tnovel\tcoding\t150\t101\t553\n")
    result = run_script("t2goAnnotationFile.py", "--classification_file", sqanti,
                        "--input_transcripts", exons, "--input_genomic", exons, "--input_protein", exons,
                        "--output", output)
    assert result.returncode == 0, result.stderr
    cds = {row[0]: row[8] for row in read_tsv(output) if row[2] == "CDS"}
    assert cds == {"PB.1.1": "ID=Protein_PB.1.1", "PB:2,1": "ID=Protein_PB_2_1"}
