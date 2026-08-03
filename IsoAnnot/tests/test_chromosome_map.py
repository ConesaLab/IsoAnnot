import os
import sys
import tempfile
import pytest

# Add scripts directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts")))

from IsoAnnot import ChromosomeMap, read_chr_ref_acc


@pytest.fixture
def sample_chr2acc_file():
    content = (
        "# Sample Chromosome to Accession Mapping Table\n"
        "1\tNC_000001.11\n"
        "2\tNC_000002.12\n"
        "X\tNC_000023.11\n"
        "scaffold_1.1\tNW_001.1\n"
        "scaffold_1.2\tNW_002.1\n"
    )
    with tempfile.NamedTemporaryFile("w+", delete=False, suffix=".tsv") as tf:
        tf.write(content)
        tf_path = tf.name
    yield tf_path
    if os.path.exists(tf_path):
        os.remove(tf_path)


def test_ensembl_leading_mode(sample_chr2acc_file):
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, leading_db="ensembl")

    # Tier 1: RefSeq -> Ensembl
    assert chr_map.get("NC_000001.11") == "1"
    assert chr_map.get("NC_000002.12") == "2"
    assert chr_map.get("NC_000023.11") == "X"

    # Tier 2: Case & Prefix normalization (Chr1, chr1 -> 1)
    assert chr_map.get("chr1") == "1"
    assert chr_map.get("Chr1") == "1"
    assert chr_map.get("CHR1") == "1"

    # Tier 3: Version-stripped RefSeq accession lookup (NC_000001 -> 1)
    assert chr_map.get("NC_000001") == "1"


def test_refseq_leading_mode(sample_chr2acc_file):
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, leading_db="refseq")

    # Tier 1: Ensembl -> RefSeq
    assert chr_map.get("1") == "NC_000001.11"
    assert chr_map.get("2") == "NC_000002.12"
    assert chr_map.get("X") == "NC_000023.11"

    # Tier 2: Normalized Ensembl name (chr1 -> NC_000001.11)
    assert chr_map.get("chr1") == "NC_000001.11"
    assert chr_map.get("Chr1") == "NC_000001.11"


def test_custom_scaffold_no_version_collision(sample_chr2acc_file):
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, leading_db="ensembl")

    # Tier 1 exact match prevents collision between scaffold_1.1 and scaffold_1.2
    assert chr_map.get("NW_001.1") == "scaffold_1.1"
    assert chr_map.get("NW_002.1") == "scaffold_1.2"


def test_unmapped_contig_fallback_and_tracking(sample_chr2acc_file):
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, leading_db="ensembl")

    # Unknown contig should return original string
    assert chr_map.get("unknown_contig_999") == "unknown_contig_999"
    assert "unknown_contig_999" in chr_map.unmapped_log
    assert chr_map.total_queries == 1


def test_len_bool_iter(sample_chr2acc_file):
    empty_map = ChromosomeMap()
    assert len(empty_map) == 0
    assert not bool(empty_map)

    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, leading_db="ensembl")
    assert len(chr_map) > 0
    assert bool(chr_map)
    assert set(chr_map) == chr_map.keys()

