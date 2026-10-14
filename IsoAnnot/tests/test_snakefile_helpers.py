import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

SNAKEFILE = os.path.join(os.path.dirname(main_path), "config", "generic", "Snakefile.smk")


def load_helpers(*names, config=None):
    """
    Copies the top-level Python functions `names` out of config/generic/Snakefile.smk and runs them in a
    namespace with `os` and `config`, so they can be tested without parsing the whole Snakefile.
    """
    lines = open(SNAKEFILE).read().splitlines()
    namespace = {"os": os, "config": config if config is not None else {}}
    for name in names:
        start = next(i for i, line in enumerate(lines) if line.startswith(f"def {name}("))
        end = start + 1
        while end < len(lines) and (not lines[end] or lines[end][0] in " \t"):
            end += 1
        exec("\n".join(lines[start:end]), namespace)
    return namespace


class TestRemoveExtension:
    """Test suite for _remove_extension."""

    @pytest.mark.parametrize("url, all, expected", [
        ("https://host/path/GCF_1_protein.faa.gz", False, "GCF_1_protein.faa"),
        ("https://host/path/Homo_sapiens.GRCh38.115.gtf.gz", True, "Homo_sapiens.GRCh38.115"),
        ("https://host/path/GCF_1_protein.faa", False, "GCF_1_protein.faa"),
        ("https://host/path/GCF_1_genomic.gtf", True, "GCF_1_genomic"),
    ])
    def test_gz_and_plain_urls(self, url, all, expected):
        """A URL without .gz keeps its name instead of crashing at parse time."""
        _remove_extension = load_helpers("_remove_extension")["_remove_extension"]
        assert _remove_extension(url, all=all) == expected
