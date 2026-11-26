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
    one per line. Returns the list of arguments isoannot.sh passed to snakemake. The stderr of the last
    run is kept in launch.stderr.
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

    def _run(*args, cwd=work_dir, **extra_env):
        result = subprocess.run(["bash", LAUNCHER, *args], cwd=cwd, env=dict(env, **extra_env),
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        _run.stderr = result.stderr
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


def test_config_drops_words_without_equals(launch):
    """Only key=value tokens after --config reach snakemake. -n or a stray word are dropped with a warning."""
    argv = launch("--database", "ensembl", "--species", "hsapiens", "--outputdir", "/abs/out",
                  "--config", "layer_go=no", "extra", "-n", "reactome=no")
    assert config_tokens(argv) == ["db=ensembl", "path_output=/abs/out", "layer_go=no", "reactome=no"]
    assert "-n" not in argv
    assert "extra" not in argv
    assert launch.stderr.splitlines() == ["Warning: --config only takes key=value pairs. Ignoring: extra",
                                          "Warning: --config only takes key=value pairs. Ignoring: -n"]


def test_config_ignores_environment(launch):
    """A CONFIG variable in the environment is not added to the --config pairs."""
    argv = launch("--database", "ensembl", "--species", "hsapiens", "--outputdir", "/abs/out",
                  "--config", "layer_go=no", CONFIG="/etc/someapp.conf")
    assert config_tokens(argv) == ["db=ensembl", "path_output=/abs/out", "layer_go=no"]


@pytest.mark.parametrize("outputdir", ["results", "./results/run1", "../results"])
def test_relative_outputdir_made_absolute(launch, tmp_path, outputdir):
    """A relative --outputdir is resolved against the current folder: run_nucimport and run_repeatmasker cd."""
    argv = launch("--database", "ensembl", "--species", "hsapiens", "--outputdir", outputdir)
    assert config_tokens(argv)[1] == f"path_output={tmp_path / 'work'}/{outputdir}"


def test_relative_outputdir_with_space_in_cwd(launch, tmp_path):
    """A space in the current folder does not split path_output into two config tokens."""
    cwd = tmp_path / "work" / "with space"
    cwd.mkdir()
    argv = launch("--database", "ensembl", "--species", "hsapiens", "--outputdir", "results", cwd=cwd)
    assert config_tokens(argv) == ["db=ensembl", f"path_output={cwd}/results"]


def test_absolute_outputdir_unchanged(launch):
    """An absolute --outputdir is passed as given."""
    argv = launch("--database", "ensembl", "--species", "hsapiens", "--outputdir", "/abs/out")
    assert config_tokens(argv)[1] == "path_output=/abs/out"


def read_version():
    from IsoAnnot import __version__
    return __version__


def test_version_option_prints_the_version_and_exits(launch):
    """--version needs no database or species, and does not start snakemake."""
    assert launch("--version") == [f"IsoAnnot {read_version()}"]


def test_version_is_the_first_line_of_a_run(launch):
    """Every run starts with the version, so the log of a job says which IsoAnnot made the output."""
    argv = launch("--database", "ensembl", "--species", "hsapiens")
    assert argv[0] == f"IsoAnnot {read_version()}"
    assert "--snakefile" in argv


def test_version_has_three_numbers():
    """isoannot.sh reads the version with a pattern, so it has to stay in the form x.y.z (or x.y.zb1)."""
    import re
    assert re.fullmatch(r"\d+\.\d+\.\d+\w*", read_version())


def test_scratch_sets_scratch_dir(launch):
    """--scratch is passed to snakemake as scratch_dir, the folder RepeatMasker works in."""
    argv = launch("--database", "ensembl", "--species", "hsapiens", "--outputdir", "/abs/out",
                  "--scratch", "/local/scratch", "--config", "layer_go=no")
    assert config_tokens(argv) == ["db=ensembl", "path_output=/abs/out", "scratch_dir=/local/scratch", "layer_go=no"]


def test_relative_scratch_made_absolute(launch, tmp_path):
    """A relative --scratch is taken from the folder isoannot.sh is run in, like --outputdir."""
    argv = launch("--database", "ensembl", "--species", "hsapiens", "--outputdir", "/abs/out", "--scratch", "scratch")
    assert config_tokens(argv) == ["db=ensembl", "path_output=/abs/out", f"scratch_dir={tmp_path}/work/scratch"]


def test_no_scratch_no_scratch_dir(launch):
    """Without --scratch no scratch_dir is passed, so the rule uses the system's temporary folder."""
    argv = launch("--database", "ensembl", "--species", "hsapiens", "--outputdir", "/abs/out")
    assert config_tokens(argv) == ["db=ensembl", "path_output=/abs/out"]
