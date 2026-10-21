import os
import sys

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from transcript2reference import group_exons_by_transcript, _retrieve_pacbio_info


def exon_line(transcript_id, start, end, strand="+", chrom="1"):
    return f'{chrom}\tPacBio\texon\t{start}\t{end}\t.\t{strand}\t.\ttranscript_id "{transcript_id}"; gene_id "PB.1";\n'


# PB.1.1 has its exons apart in the file, with those of PB.2.1 in between
GTF_LINES = [exon_line("PB.1.1", 101, 200), exon_line("PB.2.1", 1101, 1200, "-"), exon_line("PB.1.1", 301, 400),
             exon_line("PB.2.1", 1301, 1400, "-"), exon_line("PB.1.1", 501, 600)]


def test_group_exons_by_transcript():
    """Every exon of a transcript is found, also when other transcripts' exons are in between, in file order."""
    assert group_exons_by_transcript(GTF_LINES) == {
        "PB.1.1": [GTF_LINES[0], GTF_LINES[2], GTF_LINES[4]],
        "PB.2.1": [GTF_LINES[1], GTF_LINES[3]],
    }


def test_retrieve_pacbio_info():
    """
    CDS and NMD of each isoform, the same as the per-isoform scan of the whole GTF gave.
    PB.1.1: CDS 111-200 in its first exon, stop more than 55 nt before the last junction, so NMD.
    PB.2.1: minus strand, CDS from 1380 down to 1200, across the junction into the last exon, so not NMD.
    """
    orfs = {"PB.1.1\tPB.1.1.p1|30_aa|0.9|complete|11|100": "M" * 29 + "*",
            "PB.2.1\tPB.2.1.p1|27_aa|0.9|complete|21|101": "M" * 26 + "*"}
    assert _retrieve_pacbio_info(GTF_LINES, orfs) == {
        "PB.1.1": {"start": [111], "end": [200], "len": 90, "nmd": "NMD", "chr": "1"},
        "PB.2.1": {"start": [1200, 1301], "end": [1200, 1380], "len": 81, "nmd": "nonNMD", "chr": "1"},
    }
