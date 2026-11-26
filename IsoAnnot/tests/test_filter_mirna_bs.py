import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from filter_mirna_bs import mirbase_parser, site_type

MIRNA_DIR = os.path.join(os.path.dirname(main_path), "data", "global", "miRNA")

# Real duplexes from the human miRWalk 3'UTR file. Those marked "shortened" have a target structure shorter
# than the target sequence (unpaired 5' bases left out), as about 20% of the sites do.
DUPLEXES = {
    "8mer": "AGGGAGGGACGGGGGCTGTGC#GCTGGCTCCCATCTCCCTCCCA#.((((((((.((((((((.((#))))))))))...)))))))).",
    "7mer-m8": "TCACAGTGAACCGGTCTCTTT#GGGGACGGAGCACTGTGC#.(((((((..((.((((((..#))))))))..))))))).",
    "7mer-A1": "TCTGTCATTTCTTTAGGCCAATA#GCCTAGGGTGAGGGTGACAGAG#((((((((((((((((((.....#))))))))...)))))))))).",
    "6mer": "ACTGCCCTAAGTGCTCCTTCTGG#CCGGGGAGCTGAGGGGCAGG#.(((((((..(.(((((..((((#))))))))))..))))))).",
    "offset-6mer": "TGGGTCTTTGCGGGCGAGATGA#TGGTCTGGCTTTGGAAAGACCTG#.((((((((.(((((.((((.(#).)))).)))).).)))))))).",
    "seed-GU": "AGTTTTGCAGGTTTGCATCCAGC#CTGGAGGCCTCTGCAAGATTC#(((((((((((...((.(((((.#))))).)).))))))))))).",
    "noncanonical": "TGAGGTAGTAGGTTGTATGGTT#ACCTGCAGCCCCAGGCCTCC#.(((((....(((((((.(((.#))))))))))....))))).",
}
SHORTENED = {
    "8mer": "TCTGGCTCCGTGTCTTCACTCCC#GGCAGGAGGACAGGTGGAGGCGAGCCAGAT#((((((((...(((((((((.((#.))...))))))))))))))))).",
    "7mer-m8": "ATGGTTCCGTCAAGCACCATGG#CCAGCTCCAGGCAGGAAACCCCGAGAAGGAATGGTGCTTGGCAGGAACCATG#((((((((((((((((((((..#.)))))))))))).)))))))).",
    "seed-GU": "CCTGTTCTCCATTACTTGGCT#AGCTGACCAGAAGTGGCAGGAGGACAGGT#((((((((((.((((((((..#.)))..)))))..)))))))))).",
}


@pytest.mark.parametrize("expected, duplex", list(DUPLEXES.items()))
def test_site_type(expected, duplex):
    assert site_type(duplex) == expected


@pytest.mark.parametrize("expected, duplex", list(SHORTENED.items()))
def test_site_type_shortened_target_structure(expected, duplex):
    """The target structure is aligned to the 3' end of the target; without it these sites were unreadable."""
    assert site_type(duplex) == expected


@pytest.mark.parametrize("duplex", ["ACGT", "ACGTACGTA#ACGT#((((#))))", "ACGTACGTA#ACGTACGTA#((((.....#)))))))))"])
def test_site_type_unreadable(duplex):
    assert site_type(duplex) is None


@pytest.mark.parametrize("dat_file, n_mirnas", [("miRNA_high_conf.dat", 814), ("miRNA.dat", 2482)])
def test_mirbase_parser_human(dat_file, n_mirnas):
    """Human mature miRNAs with experimental evidence: the high-confidence set is the default (mirna_mirbase_set)."""
    assert len(mirbase_parser(os.path.join(MIRNA_DIR, dat_file), "Homo sapiens")) == n_mirnas


def write_inputs(tmp_path):
    """A one-miRNA miRBase file and a miRWalk table with one site of each type plus one below the score."""
    mirbase = tmp_path / "mirbase.dat"
    mirbase.write_text("ID   hsa-mir-1         standard; RNA; HSA; 80 BP.\n"
                       "DE   Homo sapiens miR-1 stem-loop\n"
                       "FT   miRNA           1..22\n"
                       "FT                   /product=\"hsa-miR-1\"\n"
                       "FT                   /evidence=experimental\n"
                       "//\n")
    rows = ["miRNA\tmRNA\tGenesymbol\tbinding_site\tduplex\tbinding_probability"]
    for n, duplex in enumerate(DUPLEXES.values()):
        rows.append(f"hsa-miR-1\tNM_{n}\tG{n}\t1,20\t{duplex}\t1.0")
    rows.append(f"hsa-miR-1\tNM_low\tGL\t1,20\t{DUPLEXES['8mer']}\t0.5")
    rows.append(f"hsa-miR-2\tNM_other\tGO\t1,20\t{DUPLEXES['8mer']}\t1.0")
    mirwalk = tmp_path / "mirwalk.txt"
    mirwalk.write_text("\n".join(rows) + "\n")
    return mirbase, mirwalk


@pytest.mark.parametrize("site_types, kept", [
    ("canonical", ["8mer", "7mer-m8", "7mer-A1"]),
    ("canonical_6mer", ["8mer", "7mer-m8", "7mer-A1", "6mer", "offset-6mer"]),
    ("all", list(DUPLEXES)),
])
def test_filter_keeps_site_types(run_script, tmp_path, site_types, kept):
    """Each --site_types setting keeps its types, with the type as a new last column; score and miRNA filters still apply."""
    mirbase, mirwalk = write_inputs(tmp_path)
    output = tmp_path / "filtered.txt"
    result = run_script("filter_mirna_bs.py", "--mirbase_file", mirbase, "--mirwalk_file", mirwalk,
                        "--species", "Homo sapiens", "--score", "0.95", "--site_types", site_types,
                        "--mirwalk_output", output)
    assert result.returncode == 0, result.stderr
    lines = output.read_text().splitlines()
    assert lines[0] == "miRNA\tmRNA\tGenesymbol\tbinding_site\tduplex\tbinding_probability\tsite_type"
    assert [line.split("\t")[6] for line in lines[1:]] == kept
    assert all(line.split("\t")[4] == DUPLEXES[line.split("\t")[6]] for line in lines[1:])
