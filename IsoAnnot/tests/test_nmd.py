import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from IsoAnnot import bases_to_last_junction, predict_nmd

# Three exons of 200, 40 and 300 nt, with introns of 800 and 1,000 nt
EXONS = [(1001, 1200), (2001, 2040), (3041, 3340)]


@pytest.mark.parametrize("cds_last, strand, expected", [
    (1150, "+", 90),      # 50 nt left in the first exon and the 40 nt of the second one
    (2010, "+", 30),      # in the second-to-last exon: the same on the genome and on the transcript
    (2040, "+", 0),       # last base before the junction
    (3100, "+", None),    # in the last exon
    (3100, "-", 99),      # minus strand: the transcript starts at 3340 and its last exon is 1001-1200
    (2010, "-", 9),
    (1100, "-", None),
])
def test_bases_to_last_junction(cds_last, strand, expected):
    """The introns are not counted: the distance is the one on the spliced transcript."""
    assert bases_to_last_junction(EXONS, cds_last, strand) == expected
    assert bases_to_last_junction(EXONS[::-1], cds_last, strand) == expected


def test_bases_to_last_junction_single_exon():
    """A transcript without junctions has no distance."""
    assert bases_to_last_junction([(1001, 1200)], 1100, "+") is None


@pytest.mark.parametrize("cds_last, strand, stop_in_cds, expected", [
    # Stop two exons before the last one, 46 nt from the junction on the transcript and 846 nt on the
    # genome, as in the main transcript of NRAS: it was always NMD with the distance on the genome
    (1194, "+", True, "nonNMD"),
    (3047, "-", True, "nonNMD"),
    # The limit: more than 50 nt after the stop codon
    (1150, "+", True, "NMD"),       # 90 nt
    (1190, "+", True, "nonNMD"),    # 50 nt
    (1189, "+", True, "NMD"),       # 51 nt
    # Ensembl and RefSeq CDS lines end before the stop codon, so its 3 bases are taken off
    (1187, "+", False, "nonNMD"),   # 53 nt after the CDS, 50 after the stop codon
    (1186, "+", False, "NMD"),      # 54 nt after the CDS, 51 after the stop codon
    (3100, "+", True, "nonNMD"),    # stop in the last exon
])
def test_predict_nmd(cds_last, strand, stop_in_cds, expected):
    """The 50 nt rule on the spliced transcript, with the same limit in the three modes."""
    assert predict_nmd(EXONS, cds_last, strand, stop_in_cds) == expected


def test_reference_sqanti_nmd_on_the_spliced_transcript(run_script, synthetic_dir, tmp_path):
    """
    referenceSQANTI.py: the CDS of TX1 ends 43 nt before the end of its first exon, and the second exon
    has 4 nt, so the stop codon is 44 nt before the last junction on the transcript and 847 nt on the
    genome. TX2 is the same with a first exon 10 nt longer: 54 nt, an NMD target.
    """
    def transcript(name, cds_end, exon_end):
        attributes = f'gene_id "G1"; transcript_id "{name}"; gene_name "G1"; transcript_biotype "protein_coding"; protein_id "P_{name}";'
        lines = [("exon", 1001, exon_end), ("exon", 2001, 2004), ("exon", 3005, 3304), ("CDS", 1101, cds_end)]
        return "".join(f"1\ttest\t{feature}\t{start}\t{end}\t.\t+\t.\t{attributes}\n" for feature, start, end in lines)

    gtf = tmp_path / "annotation.gtf"
    gtf.write_text(transcript("TX1", 1157, 1200) + transcript("TX2", 1157, 1210))
    nmd = tmp_path / "nmd.txt"
    result = run_script("referenceSQANTI.py", "--gtf_file", gtf, "--reference_file", synthetic_dir / "genome.fa",
                        "--chr_ref", synthetic_dir / "chr2acc", "--database", "ensembl",
                        "--output_classification", tmp_path / "classification.txt",
                        "--output_junctions", tmp_path / "junctions.txt", "--output_nmd", nmd)
    assert result.returncode == 0, result.stderr
    assert dict(line.split("\t") for line in nmd.read_text().splitlines()) == {"TX1": "nonNMD", "TX2": "NMD"}
