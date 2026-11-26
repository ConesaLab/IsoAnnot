import os
import subprocess
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = TESTS_DIR.parent / "scripts"
SYNTHETIC_DIR = TESTS_DIR / "fixtures" / "synthetic"


@pytest.fixture
def synthetic_dir():
    """
    Folder with the shared synthetic dataset, described in fixtures/synthetic/README.md.
    Treat it as read-only: write test outputs to tmp_path.
    """
    return SYNTHETIC_DIR


@pytest.fixture
def run_script():
    """
    Runs IsoAnnot/scripts/<name> with the test interpreter, from the scripts folder so that
    "from IsoAnnot import ..." works as in the pipeline. Returns the CompletedProcess with text
    stdout/stderr; the test checks returncode itself.
    """
    def _run(name, *args):
        return subprocess.run([sys.executable, os.path.join(SCRIPTS_DIR, name), *[str(a) for a in args]],
                              cwd=SCRIPTS_DIR, capture_output=True, text=True)
    return _run
