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


def test_ensembl_target(sample_chr2acc_file):
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, target_db="ensembl")

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


def test_refseq_target(sample_chr2acc_file):
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, target_db="refseq")

    # Tier 1: Ensembl -> RefSeq
    assert chr_map.get("1") == "NC_000001.11"
    assert chr_map.get("2") == "NC_000002.12"
    assert chr_map.get("X") == "NC_000023.11"

    # Tier 2: Normalized Ensembl name (chr1 -> NC_000001.11)
    assert chr_map.get("chr1") == "NC_000001.11"
    assert chr_map.get("Chr1") == "NC_000001.11"


def test_refseq_gtf_names_translate_to_the_ensembl_genome(sample_chr2acc_file):
    """API check for target_db="ensembl": the RefSeq chromosome bug was in the callers, see test_referenceSQANTI.py."""
    genome_names = ["1", "2", "X", "MT"]
    chr_map = read_chr_ref_acc(sample_chr2acc_file, target_db="ensembl", target_names=genome_names)
    assert [chr_map.get(c) for c in ["NC_000001.11", "NC_000002.12", "NC_000023.11"]] == ["1", "2", "X"]
    assert "NC_000001.11" in chr_map
    assert chr_map["NC_000001.11"] == "1"
    assert chr_map.unmapped_log == set()


def test_translation_missing_from_the_genome_is_unmapped(sample_chr2acc_file):
    """A name that maps but is not in the genome FASTA counts as unmapped, so the 50% abort can fire."""
    chr_map = read_chr_ref_acc(sample_chr2acc_file, target_db="ensembl", target_names=["chr1", "chr2"])
    assert chr_map.get("NC_000001.11") == "NC_000001.11"
    assert "NC_000001.11" not in chr_map
    assert chr_map.unmapped_log == {"NC_000001.11"}


def test_custom_scaffold_no_version_collision(sample_chr2acc_file):
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, target_db="ensembl")

    # Tier 1 exact match prevents collision between scaffold_1.1 and scaffold_1.2
    assert chr_map.get("NW_001.1") == "scaffold_1.1"
    assert chr_map.get("NW_002.1") == "scaffold_1.2"


def test_unmapped_contig_fallback_and_tracking(sample_chr2acc_file):
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, target_db="ensembl")

    # Unknown contig should return original string
    assert chr_map.get("unknown_contig_999") == "unknown_contig_999"
    assert "unknown_contig_999" in chr_map.unmapped_log
    assert chr_map.total_queries == 1


def test_len_bool_iter(sample_chr2acc_file):
    empty_map = ChromosomeMap()
    assert len(empty_map) == 0
    assert not bool(empty_map)

    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, target_db="ensembl")
    assert len(chr_map) > 0
    assert bool(chr_map)
    assert set(chr_map) == chr_map.keys()


def test_getitem_raises_keyerror_for_unmapped_names(sample_chr2acc_file):
    """chr_map[name] raises KeyError like a dict instead of returning the name unchanged."""
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, target_db="ensembl")
    assert chr_map["NC_000002.12"] == "2"
    with pytest.raises(KeyError):
        chr_map["unknown_contig_999"]


def test_contains_is_counted(sample_chr2acc_file):
    """"in" checks count in the failure stats; the chr_map[name] that follows them does not count twice."""
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, target_db="ensembl")
    for name in ["NC_000001.11", "MT", "MT"]:
        if name in chr_map:
            chr_map[name]
    assert (chr_map.total_queries, chr_map.unmapped_queries, chr_map.unmapped_log) == (3, 2, {"MT"})


def test_check_unmapped_counts_lookups_not_names(sample_chr2acc_file):
    """Three lookups of one unmapped name out of four is 75%, not 1/4 distinct names: the abort fires."""
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, target_db="ensembl")
    for name in ["NC_000001.11", "MT", "MT", "MT"]:
        chr_map.get(name, name)
    with pytest.raises(RuntimeError, match=r"test.py: 3/4 chromosome lookups failed .*unmapped names: MT\)"):
        chr_map.check_unmapped("test.py")


def test_check_unmapped_accepts_a_few_unmapped_names(sample_chr2acc_file):
    """MT and scaffolds missing from chr2acc are normal: half or less of the lookups does not stop the run."""
    chr_map = ChromosomeMap(mapping_file=sample_chr2acc_file, target_db="ensembl")
    for name in ["NC_000001.11", "MT"]:
        chr_map.get(name, name)
    chr_map.check_unmapped("test.py")
    assert (chr_map.total_queries, chr_map.unmapped_queries) == (2, 1)

