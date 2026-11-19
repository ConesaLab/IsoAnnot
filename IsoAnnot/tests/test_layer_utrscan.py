import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from layer_utrscan import UTRSITE_FILE, get_motif_feature, read_utrsite_regions


def test_read_utrsite_regions():
    """The catalogue IsoAnnot ships gives the UTR of each of its 59 motifs."""
    regions = read_utrsite_regions(UTRSITE_FILE)
    assert len(regions) == 59
    assert {motif: regions[motif] for motif in ("uORF", "IRES", "TOP", "PAS", "MBE", "K-BOX", "IRE", "SXL_BS")} == {
        "uORF": "5", "IRES": "5", "TOP": "5", "PAS": "3", "MBE": "3", "K-BOX": "3", "IRE": "all", "SXL_BS": "all"}
    assert sorted(set(regions.values())) == ["3", "5", "all"]


@pytest.mark.parametrize("motif, region, start, end, expected", [
    ("uORF", "5", 20, 90, "uORF"),          # before the CDS (101-550)
    ("IRES", "5", 20, 90, "5UTRmotif"),
    ("IRE", "all", 20, 90, "5UTRmotif"),
    ("MBE", "3", 20, 26, None),             # a 3'UTR motif found in the 5'UTR
    ("uORF", "5", 600, 700, None),          # an "upstream ORF" after the CDS
    ("IRES", "5", 600, 700, None),
    ("MBE", "3", 600, 606, "3UTRmotif"),
    ("IRE", "all", 600, 630, "3UTRmotif"),
    ("PAS", "3", 720, 750, "PAS"),
    ("uORF", "5", 90, 160, None),           # across the start of the CDS
    ("MBE", "3", 300, 306, None),           # inside the CDS
])
def test_get_motif_feature_coding(motif, region, start, end, expected):
    """In a coding transcript a motif is kept only in the UTR of its catalogue entry."""
    assert get_motif_feature(motif, region, start, end, 101, 550, keep_noncoding=True) == expected


@pytest.mark.parametrize("motif, region, keep_noncoding, expected", [
    ("MBE", "3", True, "ncRNAmotif"),
    ("IRE", "all", True, "ncRNAmotif"),
    ("PAS", "3", True, "PAS"),
    ("uORF", "5", True, None),      # no CDS to be upstream of
    ("IRES", "5", True, None),
    ("MBE", "3", False, None),
    ("PAS", "3", False, None),
])
def test_get_motif_feature_noncoding(motif, region, keep_noncoding, expected):
    """A transcript without a CDS has no UTR: its motifs are ncRNAmotif, or left out when asked."""
    assert get_motif_feature(motif, region, 300, 330, None, None, keep_noncoding=keep_noncoding) == expected


CLASSIFICATION = ("isoform\tchrom\tstrand\tlength\tstructural_category\tassociated_gene\tassociated_transcript\tRTS_stage"
                  "\tcoding\tORF_length\tCDS_start\tCDS_end\n"
                  "TX1\t1\t+\t750\tprotein_coding\tSYNA\tTX1\tNA\tcoding\t150\t101\t550\n"
                  "NC1\t1\t-\t700\tlncRNA\tSYNC\tNC1\tNA\tnon_coding\tNA\tNA\tNA\n")
UTRSCAN = ("U0033 - uORF            - Upstream Open Reading Frame (uORF)\n"
           "U0035 - MBE             - Musashi binding element (MBE)\n"
           "U0043 - PAS             - Polyadenylation Signal (PAS)\n"
           "U9999 - NEWMOTIF        - A motif that is not in the catalogue (NEWMOTIF)\n"
           "TX1 length=750 : uORF [20,90]\n"
           "TX1 length=750 : MBE [30,36]\n"
           "TX1 length=750 : uORF [600,700]\n"
           "TX1 length=750 : MBE [610,616]\n"
           "TX1 length=750 : PAS [720,750]\n"
           "TX1 length=750 : NEWMOTIF [40,50]\n"
           "NC1 length=700 : uORF [100,200]\n"
           "NC1 length=700 : MBE [300,306]\n"
           "NC1 length=700 : PAS [670,700]\n")


def run_layer(run_script, tmp_path, *options):
    classification = tmp_path / "classification.txt"
    classification.write_text(CLASSIFICATION)
    utrscan = tmp_path / "utrscan.txt"
    utrscan.write_text(UTRSCAN)
    output = tmp_path / "layer_utrscan.gtf"
    result = run_script("layer_utrscan.py", "--utrscan_file", utrscan, "--classification_file", classification,
                        "--output", output, *options)
    assert result.returncode == 0, result.stderr
    return [tuple(line.split("\t")[:5]) + (line.split("\t")[8].split(";")[1].strip(),) for line in output.read_text().splitlines()], result.stderr


def test_layer_utrscan_keeps_each_motif_in_its_utr(run_script, tmp_path):
    """
    UtrScan reports every motif along the whole transcript. The upstream ORF after the CDS and the 3'UTR
    motif before it are left out, and the motifs of the non-coding transcript, which used to be 3UTRmotif,
    are ncRNAmotif. A motif that is not in the catalogue is kept where it is found.
    """
    rows, log = run_layer(run_script, tmp_path)
    assert rows == [
        ("TX1", "UTRsite", "uORF", "20", "90", "Name=uORF"),
        ("TX1", "UTRsite", "3UTRmotif", "610", "616", "Name=MBE"),
        ("TX1", "UTRsite", "PAS", "720", "750", "Name=PAS"),
        ("TX1", "UTRsite", "5UTRmotif", "40", "50", "Name=NEWMOTIF"),
        ("NC1", "UTRsite", "ncRNAmotif", "300", "306", "Name=MBE"),
        ("NC1", "UTRsite", "PAS", "670", "700", "Name=PAS"),
    ]
    assert ("UtrScan hits: 1 3UTRmotif, 1 5UTRmotif, 2 PAS, 2 left out, coding, 1 left out, non-coding, "
            "1 ncRNAmotif, 1 uORF") in log


def test_layer_utrscan_skip_noncoding(run_script, tmp_path):
    """--skip_noncoding (utrscan_noncoding: False) leaves out every motif of the transcripts without a CDS."""
    rows, _ = run_layer(run_script, tmp_path, "--skip_noncoding")
    assert [row[0] for row in rows] == ["TX1", "TX1", "TX1", "TX1"]


def test_rename_features_ncrna_motif(run_script, tmp_path):
    """ncRNAmotif is a positional transcript feature, as the other UTRsite ones."""
    gff3 = tmp_path / "annotation.gff3"
    gff3.write_text("NC1\tUTRsite\tncRNAmotif\t300\t306\t.\t-\t.\tID=U0035; Name=MBE; Desc=Musashi binding element (MBE)\n")
    result = run_script("renameFeatures.py", gff3)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "annotation.gff3_mod").read_text() == (
        "NC1\tUTRsite\tncRNAmotif\t300\t306\t.\t-\t.\tID=U0035; Name=MBE; Desc=Musashi binding element (MBE); PosType=T\n")
    assert "does not know" not in result.stderr
