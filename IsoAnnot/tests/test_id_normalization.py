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


CONFIG_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "config"))

# (isoform, refProt, transcript ID in the tool output, protein ID in the tool output)
ID_CASES = {
    "ensembl": ("ENST00000000001", "ENSP00000000001", "ENST00000000001.4", "ENSP00000000001.2"),
    "refseq": ("NM_900001.1", "NP_900001.1", "NM_900001.1", "NP_900001.1"),
    "mytranscripts": ("PB.1.1", None, "PB.1.1", "PB.1.1"),
}


def write_classification(path, isoform, ref_prot):
    """One coding transcript, with a refProt column only when there is a reference protein (as in SQANTI3)."""
    header = ["isoform", "chrom", "strand", "length", "structural_category", "associated_gene",
              "associated_transcript", "RTS_stage", "coding", "ORF_length", "CDS_start", "CDS_end"]
    values = [isoform, "1", "+", "750", "full-splice_match", "SYNA", "ENST00000000001", "FALSE", "coding",
              "150", "101", "550"]
    if ref_prot:
        header.append("refProt")
        values.append(ref_prot)
    path.write_text("\t".join(header) + "\n" + "\t".join(values) + "\n")


@pytest.mark.parametrize("keep_version", ["False", "True"])
@pytest.mark.parametrize("mode", ID_CASES)
def test_layer_repeatmasker_ids(tmp_path, run_script, mode, keep_version):
    """RefSeq and PacBio IDs were cut at the first dot with keep_version False, and Ensembl IDs never matched with True."""
    isoform, ref_prot, transcript_id, _ = ID_CASES[mode]
    classification = tmp_path / "classification.txt"
    write_classification(classification, isoform, ref_prot)
    repeatmasker_file = tmp_path / "transcripts.fa.out"
    repeatmasker_file.write_text(
        "   SW   perc perc perc  query      position in query    matching  repeat\n"
        "score   div. del. ins.  sequence   begin  end  (left)   repeat    class/family\n"
        "\n"
        f"  463   1.3  0.6  1.7  {transcript_id}    1   155 (595) +  AluY  SINE/Alu  1  155 (0)  1\n")
    output = tmp_path / "layer_repeatmasker.gtf"
    result = run_script("layer_repeatmasker.py", "--keep_version", keep_version, "--repeatmasker_file", repeatmasker_file,
                        "--classification_file", classification, "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [
        f"{isoform}\tRepeatMasker\trepeat\t1\t155\t.\t+\t.\tID=SINE/Alu; Name=SINE/Alu; Desc=AluY"]


@pytest.mark.parametrize("keep_version", ["False", "True"])
@pytest.mark.parametrize("mode", ID_CASES)
def test_layer_utrscan_ids(tmp_path, run_script, mode, keep_version):
    """RefSeq and PacBio IDs were cut at the first dot with keep_version False, and Ensembl IDs never matched with True."""
    isoform, ref_prot, transcript_id, _ = ID_CASES[mode]
    classification = tmp_path / "classification.txt"
    write_classification(classification, isoform, ref_prot)
    utrscan_file = tmp_path / "utrscan.txt"
    utrscan_file.write_text(
        "U0001 - uORF - Upstream Open Reading Frame\n"
        f"{transcript_id} length=750 : uORF [20,40]\n")
    output = tmp_path / "layer_utrscan.gtf"
    result = run_script("layer_utrscan.py", "--keep_version", keep_version, "--utrscan_file", utrscan_file,
                        "--classification_file", classification, "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [
        f"{isoform}\tUTRsite\tuORF\t20\t40\t.\t+\t.\tID=U0001; Name=uORF; Desc=Upstream Open Reading Frame"]


@pytest.mark.parametrize("keep_version", ["False", "True"])
@pytest.mark.parametrize("mode", ID_CASES)
def test_layer_interproscan_ids(tmp_path, run_script, mode, keep_version):
    """RefSeq and PacBio protein IDs were cut at the first dot with keep_version False, and Ensembl IDs never matched with True."""
    isoform, ref_prot, _, protein_id = ID_CASES[mode]
    classification = tmp_path / "classification.txt"
    write_classification(classification, isoform, ref_prot)
    interproscan_file = tmp_path / "interproscan.tsv"
    interproscan_file.write_text(f"{protein_id}\tPfam\tx\tPF00001\t7tm_1\t7 transmembrane receptor\t10\t50\n")
    pfam_file = tmp_path / "pfam.tsv"
    pfam_file.write_text("")
    output = tmp_path / "layer_interproscan.gtf"
    result = run_script("layer_interproscan.py", "--keep_version", keep_version, "--interproscan_file", interproscan_file,
                        "--pfam_file", pfam_file, "--classification_file", classification, "--output", output)
    assert result.returncode == 0, result.stderr
    assert output.read_text().splitlines() == [
        f"{isoform}\tPFAM\tDOMAIN\t10\t50\t.\t.\t.\tID=PF00001; Name=7tm_1; Desc=7 transmembrane receptor"]


def test_get_fasta_sequences_keeps_pacbio_ids(tmp_path):
    """With keep_version False, PB.1.1 and PB.1.2 both became the key PB and one sequence was lost."""
    pytest.importorskip("Bio")
    from uniprotPhosphosite_genomicCoordinates import get_fasta_sequences
    fasta = tmp_path / "orf.faa"
    fasta.write_text(">PB.1.1 gene=x\nMKRK\n>PB.1.2 gene=x\nMAAA\n>ENSP00000000001.2 pep\nMCCC\n")
    # Same regex as uniprotPhosphosite_annotation.py uses for the ORF FASTA.
    assert get_fasta_sequences(fasta_files=[str(fasta)], matching_regex="(.*?)\\.?\\s", keep_version="False") == {
        "PB.1.1": "MKRK", "PB.1.2": "MAAA", "ENSP00000000001": "MCCC"}


@pytest.mark.parametrize("config_name", [
    f"{db}/{species}"
    for db in ("refseq", "mytranscripts")
    for species in sorted(os.listdir(os.path.join(CONFIG_DIR, db)))
    if os.path.isfile(os.path.join(CONFIG_DIR, db, species, "config.yaml"))
])
def test_versioned_configs(config_name):
    """RefSeq and mytranscripts IDs are always versioned, so transcript_versioned must be True."""
    with open(os.path.join(CONFIG_DIR, config_name, "config.yaml")) as handle:
        values = [line.split(":", 1)[1].strip() for line in handle if line.startswith("transcript_versioned:")]
    assert values == ["True"]
