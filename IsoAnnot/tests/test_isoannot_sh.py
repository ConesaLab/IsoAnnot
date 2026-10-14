import os
import subprocess
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

ISOANNOT_DIR = os.path.dirname(main_path)
LAUNCHER = os.path.join(ISOANNOT_DIR, "isoannot.sh")


@pytest.fixture
def launch(tmp_path):
    """
    Runs isoannot.sh from tmp_path/work with a fake "snakemake" first on PATH that prints its arguments,
    one per line. Returns the list of arguments isoannot.sh passed to snakemake.
    """
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake = bin_dir / "snakemake"
    fake.write_text('#!/bin/sh\nfor a in "$@"; do printf "%s\\n" "$a"; done\n')
    fake.chmod(0o755)
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    env = dict(os.environ, PATH=f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    env.pop("SLURM_CPUS_PER_TASK", None)

    def _run(*args, **extra_env):
        result = subprocess.run(["bash", LAUNCHER, *args], cwd=work_dir, env=dict(env, **extra_env),
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        return result.stdout.splitlines()
    return _run


def config_tokens(argv):
    """The key=value tokens isoannot.sh passed after --config, up to the next option."""
    start = argv.index("--config") + 1
    end = argv.index("--directory")
    return argv[start:end]


@pytest.mark.parametrize("config_args", [
    ["--config", "layer_go=no", "reactome=no"],
    ["--config", "layer_go=no reactome=no"],
    ["--config", "layer_go=no", "--config", "reactome=no"],
])
def test_config_keeps_every_pair(launch, config_args):
    """--config takes several key=value pairs, unquoted, quoted as one string or repeated."""
    argv = launch("--database", "ensembl", "--species", "hsapiens", "--outputdir", "/abs/out", *config_args)
    assert config_tokens(argv) == ["db=ensembl", "path_output=/abs/out", "layer_go=no", "reactome=no"]


def test_config_stops_at_next_option(launch):
    """Options after the --config pairs are still read by isoannot.sh, not passed as config."""
    argv = launch("--config", "layer_go=no", "reactome=no", "--database", "ensembl", "--species", "hsapiens",
                  "--outputdir", "/abs/out")
    assert config_tokens(argv) == ["db=ensembl", "path_output=/abs/out", "layer_go=no", "reactome=no"]
    assert argv[argv.index("--snakefile") + 1] == os.path.join(ISOANNOT_DIR, "config", "ensembl", "hsapiens", "Snakefile.smk")


def test_config_ignores_environment(launch):
    """A CONFIG variable in the environment is not added to the --config pairs."""
    argv = launch("--database", "ensembl", "--species", "hsapiens", "--outputdir", "/abs/out",
                  "--config", "layer_go=no", CONFIG="/etc/someapp.conf")
    assert config_tokens(argv) == ["db=ensembl", "path_output=/abs/out", "layer_go=no"]
