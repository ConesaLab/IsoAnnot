import os
import sys

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from layer_interproscan import get_clan_descriptions

PFAM_CLANS = ("PF00028\tCL0159\tE-set\tCadherin\tCadherin domain\n"
              "PF01825\t\t\tGPS\tGPCR proteolysis site, GPS, motif\n")
PFAM_C = ("# STOCKHOLM 1.0\n#=GF ID   E-set\n#=GF AC   CL0159.20\n#=GF DE   Ig-like fold superfamily (E-set)\n"
          "#=GF AU   Someone\n//\n#=GF ID   LRR\n#=GF AC   CL0022.41\n#=GF DE   Leucine Rich Repeat\n//\n")


def run_layer(run_script, synthetic_dir, tmp_path, interproscan_rows, with_clan_file=True):
    interproscan = tmp_path / "interproResultsPfam.tsv"
    interproscan.write_text(interproscan_rows)
    pfam = tmp_path / "Pfam-A.clans.tsv"
    pfam.write_text(PFAM_CLANS)
    clans = tmp_path / "Pfam-C"
    clans.write_text(PFAM_C)
    output = tmp_path / "layer_interproscan.gtf"
    args = ["--interproscan_file", interproscan, "--pfam_file", pfam, "--classification_file", synthetic_dir / "ensembl_sqanti_classification.txt",
            "--output", output]
    if with_clan_file:
        args += ["--pfam_clan_desc", clans]
    result = run_script("layer_interproscan.py", *args)
    assert result.returncode == 0, result.stderr
    return output.read_text().splitlines()


def test_get_clan_descriptions():
    """Clan accessions lose their version, so they match the Pfam-A.clans.tsv column."""
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as h:
        h.write(PFAM_C)
    assert get_clan_descriptions(h.name) == {"CL0159": "Ig-like fold superfamily (E-set)", "CL0022": "Leucine Rich Repeat"}
    os.remove(h.name)


def test_clan_line_names_the_clan(run_script, synthetic_dir, tmp_path):
    """The clan line has the clan's name and description, as in the legacy tappAS files, not the family's."""
    lines = run_layer(run_script, synthetic_dir, tmp_path, "ENSP00000000001\tPfam\t-\tPF00028\tCadherin\tCadherin domain\t250\t343\n")
    assert lines == [
        "ENST00000000001\tPFAM\tDOMAIN\t250\t343\t.\t.\t.\tID=PF00028; Name=Cadherin; Desc=Cadherin domain",
        "ENST00000000001\tPFAM\tclan\t.\t.\t.\t.\t.\tID=CL0159; Name=E-set; Desc=Ig-like fold superfamily (E-set)",
    ]


def test_no_clan_line_for_a_family_without_clan(run_script, synthetic_dir, tmp_path):
    """A family with no clan gets its domain line only, not a clan line with an empty ID."""
    lines = run_layer(run_script, synthetic_dir, tmp_path, "ENSP00000000001\tPfam\t-\tPF01825\tGPS\tGPCR proteolysis site, GPS, motif\t2409\t2454\n")
    assert lines == ["ENST00000000001\tPFAM\tDOMAIN\t2409\t2454\t.\t.\t.\tID=PF01825; Name=GPS; Desc=GPCR proteolysis site, GPS, motif"]


def test_clan_description_without_clan_file(run_script, synthetic_dir, tmp_path):
    """Without Pfam's clan file the clan name is used as description, instead of the family's."""
    lines = run_layer(run_script, synthetic_dir, tmp_path, "ENSP00000000001\tPfam\t-\tPF00028\tCadherin\tCadherin domain\t250\t343\n", with_clan_file=False)
    assert lines[1] == "ENST00000000001\tPFAM\tclan\t.\t.\t.\t.\t.\tID=CL0159; Name=E-set; Desc=E-set"


def test_default_description_by_database(run_script, synthetic_dir, tmp_path):
    """An InterProScan hit without description gets its database's default text; it was always None."""
    lines = run_layer(run_script, synthetic_dir, tmp_path, "ENSP00000000001\tMOBIDB_LITE\tDOMAIN\tmobidb-lite\tconsensus disorder\tNone\t1\t20\n")
    assert lines == ["ENST00000000001\tMOBIDB_LITE\tDISORDER\t1\t20\t.\t.\t.\tID=mobidb-lite; Name=consensus disorder; Desc=consensus disorder prediction"]


def test_get_clan_descriptions_with_stray_bytes(tmp_path):
    """The real Pfam-C has a non-UTF-8 byte (0x96) in a comment line; it must not stop the layer."""
    import gzip
    clan_file = tmp_path / "Pfam-C.gz"
    with gzip.open(clan_file, "wb") as h:
        h.write(b"#=GF ID   LRR\n#=GF AC   CL0022.41\n#=GF DE   Leucine Rich Repeat\n#=GF CC   cytochrome c\x96dependent\n//\n")
    assert get_clan_descriptions(str(clan_file)) == {"CL0022": "Leucine Rich Repeat"}
