import json
import os
import subprocess
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

SNAKEFILE = os.path.join(os.path.dirname(main_path), "config", "generic", "Snakefile.smk")
SNAKEMAKE = os.environ.get("ISOANNOT_SNAKEMAKE", "/media/pabloati/HDD1/conda/envs/snakemake/bin/snakemake")


def load_helpers(*names, config=None, **variables):
    """
    Copies the top-level Python functions `names` out of config/generic/Snakefile.smk and runs them in a
    namespace with `os`, `config` and the extra `variables` (e.g. prefix, path_output), so they can be
    tested without parsing the whole Snakefile.
    """
    namespace = {"os": os, "config": config if config is not None else {}, **variables}
    for name in names:
        exec(helper_source(name), namespace)
    return namespace


def helper_source(name):
    """Source of the top-level function `name` in config/generic/Snakefile.smk."""
    lines = open(SNAKEFILE).read().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(f"def {name}("))
    end = start + 1
    while end < len(lines) and (not lines[end] or lines[end][0] in " \t"):
        end += 1
    return "\n".join(lines[start:end])


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


class TestGatherInterproscan:
    """Test suite for gather_interproscan."""

    # Runs in the Python of the Snakemake env, for the real expand and glob_wildcards. The checkpoint
    # object is replaced by one whose output folder is xml_dir.
    RUNNER = """
import json, os, sys
from types import SimpleNamespace
from snakemake.io import expand, glob_wildcards
xml_dir = sys.argv[1]
checkpoints = SimpleNamespace(run_interproscan=SimpleNamespace(get=lambda **wildcards: SimpleNamespace(output=[xml_dir])))
namespace = {"os": os, "expand": expand, "glob_wildcards": glob_wildcards, "checkpoints": checkpoints}
exec(sys.stdin.read(), namespace)
print(json.dumps(namespace["gather_interproscan"]({"db": "ensembl"})))
"""

    def run_gather(self, xml_dir):
        python = os.path.join(os.path.dirname(SNAKEMAKE), "python")
        if not os.path.isfile(python):
            pytest.skip(f"Python of the Snakemake env not found at {python} (set ISOANNOT_SNAKEMAKE to use another one)")
        return subprocess.run([python, "-c", self.RUNNER, str(xml_dir)], input=helper_source("gather_interproscan"),
                              capture_output=True, text=True)

    def gather(self, xml_dir):
        result = self.run_gather(xml_dir)
        assert result.returncode == 0, result.stderr
        return json.loads(result.stdout)

    def test_all_xml_files(self, tmp_path):
        """Every InterProScan XML is passed to parse_interproscan, not only the first one."""
        for name in ["proteins_1.xml", "proteins_2.xml", "notes.txt"]:
            (tmp_path / name).write_text("")
        assert sorted(self.gather(tmp_path)) == [str(tmp_path / "proteins_1.xml"), str(tmp_path / "proteins_2.xml")]

    def test_no_xml_files(self, tmp_path):
        """No XML gives an error that names the InterProScan folder, instead of an IndexError."""
        result = self.run_gather(tmp_path)
        assert result.returncode == 1
        assert f"ValueError: InterProScan wrote no XML file in {tmp_path}. Please check run_interproscan.log in the logs folder." in result.stderr
