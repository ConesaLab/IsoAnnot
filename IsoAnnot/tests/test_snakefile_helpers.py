import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

SNAKEFILE = os.path.join(os.path.dirname(main_path), "config", "generic", "Snakefile.smk")


def load_helpers(*names, config=None, **variables):
    """
    Copies the top-level Python functions `names` out of config/generic/Snakefile.smk and runs them in a
    namespace with `os`, `config` and the extra `variables` (e.g. prefix, path_output), so they can be
    tested without parsing the whole Snakefile.
    """
    lines = open(SNAKEFILE).read().splitlines()
    namespace = {"os": os, "config": config if config is not None else {}, **variables}
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


class TestOutputLayerDb:
    """Test suite for _output_layer_db."""

    @pytest.mark.parametrize("layer_go, expected", [
        ("si", "/out/data/Hsapiens/output/{db}/layers/layer_go.gtf"),
        (True, "/out/data/Hsapiens/output/{db}/layers/layer_go.gtf"),
        ("yes", "/out/data/Hsapiens/output/{db}/layers/layer_go.gtf"),
        ("no", []),
        (False, []),
        ("", []),
    ])
    def test_layer_go_gate(self, layer_go, expected):
        """Without an external rule, the layer_go value alone turns the layer on or off; "no" no longer counts as on."""
        helpers = load_helpers("is_feature_enabled", "_output_layer_db", config={"layer_go": layer_go},
                               prefix="Hsapiens", db="ensembl", path_output="/out")
        assert helpers["_output_layer_db"]("layer_go", []) == expected

    @pytest.mark.parametrize("layer_go", ["si", True, "yes", "no", False])
    def test_layer_go_values(self, layer_go):
        """A boolean layer_go (YAML true/false) no longer raises TypeError at parse time."""
        helpers = load_helpers("is_feature_enabled", "_output_layer_db", config={"layer_go": layer_go},
                               prefix="Hsapiens", db="ensembl", path_output="/out")
        assert helpers["_output_layer_db"]("layer_go", lambda x: "hsapiens_gene_ensembl") == "/out/data/Hsapiens/output/{db}/layers/layer_go.gtf"
