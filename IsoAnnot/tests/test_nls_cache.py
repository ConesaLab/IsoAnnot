import os
import subprocess
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from nls_cache import NUCIMPORT_HEADER, NUCIMPORT_STARS, STORE_COLUMNS, sequence_hash

SCRIPT = os.path.join(main_path, "nls_cache.py")

SEQ_A = "MSPKKKRKVEDPAAAALLLLKRKRRSPPKKKRKVAAAGGGLLL"
SEQ_B = "MKRPRGRPKKKRKVLLAAAPEDKRRKAAAGGGLLLAAAGGGLLLDDE"
SEQ_C = "MAAAAGGGGLLLLVVVVSSSSAAAAGGGGLLLLVVVVSSSS"
VALUES_A = ["0.9949979710297792", "0.806466703423478", "AALLLLKRKRRSPPKKKRKV", "Class1", "15"]
VALUES_B = ["0.9933396900000001", "0.70", "PRGRPKKKRKVLLAAAPEDK", "Class1", "4"]
VALUES_C = ["0.0031678638823667", "9.177442751606701E-5", "AAAAGGGGLLLLVVVVSSSS", "Class2", "2"]


def nucimport_line(prot_id, values):
    """One NucImport result line: the jar ends every line with a tab."""
    return "\t".join([prot_id] + values) + "\t\n"


def nucimport_output(*rows):
    return f"{NUCIMPORT_HEADER}\n{NUCIMPORT_STARS}\n" + "".join(nucimport_line(p, v) for p, v in rows)


@pytest.fixture
def jar(tmp_path):
    jar = tmp_path / "NucImport.jar"
    jar.write_bytes(b"jar v1")
    return jar


@pytest.fixture
def cache(tmp_path, jar):
    """Runs nls_cache.py lookup or store with the cache in tmp_path/cache. Returns the CompletedProcess."""
    def _run(command, *args, model="Model4", jar_path=None, cache_option=None):
        cmd = [sys.executable, SCRIPT, command, *[str(a) for a in args],
               "--model", model, "--organism", "Mouse", "--id_flag", "ID=F", "--jar", str(jar_path or jar)]
        cmd += cache_option or ["--cache_dir", str(tmp_path / "cache")]
        return subprocess.run(cmd, capture_output=True, text=True)
    return _run


def lookup(cache, tmp_path, fasta_text, **kwargs):
    """Lookup of the proteins in fasta_text. Returns the CompletedProcess, the to_predict text and the cached text."""
    fasta = tmp_path / "unique_proteins.fa"
    fasta.write_text(fasta_text)
    to_predict, cached = tmp_path / "to_predict.fa", tmp_path / "cached_results.txt"
    result = cache("lookup", "--fasta", fasta, "--to_predict", to_predict, "--cached", cached, **kwargs)
    assert result.returncode == 0, result.stderr
    return result, to_predict.read_text(), cached.read_text()


def store(cache, tmp_path, fasta_text, *outputs, **kwargs):
    """Store of the NucImport outputs for the proteins in fasta_text. Returns the CompletedProcess."""
    fasta = tmp_path / "sent.fa"
    fasta.write_text(fasta_text)
    results = []
    for n, text in enumerate(outputs):
        results.append(tmp_path / f"output_chunk_{n}.txt")
        results[-1].write_text(text)
    return cache("store", "--fasta", fasta, "--results", *results, "--name", "Hsapiens_ensembl",
                 "--missing", tmp_path / "missing.txt", **kwargs)


def store_files(tmp_path):
    cache_dir = tmp_path / "cache"
    return sorted(os.listdir(cache_dir)) if cache_dir.exists() else []


def test_lookup_empty_cache(cache, tmp_path):
    """Without a store every protein is predicted; the cached file has only the NucImport header."""
    result, to_predict, cached = lookup(cache, tmp_path, f">P1\n{SEQ_A}\n>P2\n{SEQ_B}\n")
    assert to_predict == f">P1\n{SEQ_A}\n>P2\n{SEQ_B}\n"
    assert cached == nucimport_output()
    assert "0 of 2 proteins found in the cache, 2 to predict" in result.stderr


def test_store_then_lookup_other_id(cache, tmp_path):
    """A sequence stored under one ID is found under another one, e.g. another representative or database."""
    result = store(cache, tmp_path, f">ENSP1.1\n{SEQ_A}\n>ENSP2.1\n{SEQ_B}\n",
                   nucimport_output(("ENSP1.1", VALUES_A)), nucimport_output(("ENSP2.1", VALUES_B)))
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "missing.txt").read_text() == ""
    [store_file] = store_files(tmp_path)
    assert store_file.startswith("Hsapiens_ensembl_")
    lines = (tmp_path / "cache" / store_file).read_text().splitlines()
    assert lines[0] == "\t".join(STORE_COLUMNS)
    assert lines[1].split("\t")[:3] == [sequence_hash(SEQ_A), "Model4", "Mouse"]

    _, to_predict, cached = lookup(cache, tmp_path, f">NP_9.1\n{SEQ_B}\n>PB.3.1\n{SEQ_C}\n>NP_7.2\n{SEQ_A}\n")
    assert to_predict == f">PB.3.1\n{SEQ_C}\n"
    assert cached == nucimport_output(("NP_9.1", VALUES_B), ("NP_7.2", VALUES_A))


@pytest.mark.parametrize("option", ["model", "jar"])
def test_lookup_other_settings_is_miss(cache, tmp_path, option):
    """A result is only reused with the same model and the same jar."""
    store(cache, tmp_path, f">P1\n{SEQ_A}\n", nucimport_output(("P1", VALUES_A)))
    if option == "model":
        kwargs = {"model": "Model2"}
    else:
        new_jar = tmp_path / "NucImport_v2.jar"
        new_jar.write_bytes(b"jar v2")
        kwargs = {"jar_path": new_jar}
    _, to_predict, cached = lookup(cache, tmp_path, f">P1\n{SEQ_A}\n", **kwargs)
    assert to_predict == f">P1\n{SEQ_A}\n"
    assert cached == nucimport_output()


def test_store_only_new_results(cache, tmp_path):
    """Results already in the store are not written again, and nothing new writes no file."""
    store(cache, tmp_path, f">P1\n{SEQ_A}\n", nucimport_output(("P1", VALUES_A)))
    result = store(cache, tmp_path, f">P1\n{SEQ_A}\n", nucimport_output(("P1", VALUES_A)))
    assert result.returncode == 0, result.stderr
    assert "No new NucImport results to store" in result.stderr
    assert len(store_files(tmp_path)) == 1


def test_store_some_missing(cache, tmp_path):
    """A protein without a result is listed and left out of the store, so the next run predicts it again."""
    result = store(cache, tmp_path, f">P1\n{SEQ_A}\n>P2\n{SEQ_B}\n>P3\n{SEQ_C}\n",
                   nucimport_output(("P1", VALUES_A), ("P3", VALUES_C)))
    assert result.returncode == 0, result.stderr
    assert "WARNING: 1 of 3 proteins sent to NucImport got no result" in result.stderr
    assert (tmp_path / "missing.txt").read_text() == "P2\n"
    _, to_predict, _ = lookup(cache, tmp_path, f">P1\n{SEQ_A}\n>P2\n{SEQ_B}\n>P3\n{SEQ_C}\n")
    assert to_predict == f">P2\n{SEQ_B}\n"


def test_store_most_missing_stops(cache, tmp_path):
    """More than half of the proteins without a result means NucImport did not work: the run stops and nothing is stored."""
    result = store(cache, tmp_path, f">P1\n{SEQ_A}\n>P2\n{SEQ_B}\n>P3\n{SEQ_C}\n",
                   nucimport_output(("P1", VALUES_A)), "")
    assert result.returncode == 1
    assert "ERROR: 2 of 3 proteins sent to NucImport got no result" in result.stderr
    assert (tmp_path / "missing.txt").read_text() == "P2\nP3\n"
    assert store_files(tmp_path) == []


def test_store_nothing_sent(cache, tmp_path):
    """Every protein cached: nothing was sent to NucImport, every chunk is empty and the run goes on."""
    result = store(cache, tmp_path, "", "", "")
    assert result.returncode == 0, result.stderr
    assert store_files(tmp_path) == []


def test_no_cache(cache, tmp_path):
    """With the cache off, lookup sends every protein and store checks the results without storing them."""
    store(cache, tmp_path, f">P1\n{SEQ_A}\n", nucimport_output(("P1", VALUES_A)))
    _, to_predict, cached = lookup(cache, tmp_path, f">P1\n{SEQ_A}\n", cache_option=["--no_cache"])
    assert to_predict == f">P1\n{SEQ_A}\n"
    assert cached == nucimport_output()
    result = store(cache, tmp_path, f">P2\n{SEQ_B}\n", "", cache_option=["--no_cache"])
    assert result.returncode == 1
    assert len(store_files(tmp_path)) == 1


def test_store_conflict_first_file_wins(cache, tmp_path):
    """Two store files that disagree on a sequence: the first in name order is used, with a warning."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    jar_md5 = __import__("hashlib").md5(b"jar v1").hexdigest()
    for name, values in [("a.tsv", VALUES_A), ("b.tsv", VALUES_C)]:
        row = [sequence_hash(SEQ_A), "Model4", "Mouse", "ID=F", jar_md5] + values
        (cache_dir / name).write_text("\t".join(STORE_COLUMNS) + "\n" + "\t".join(row) + "\n")
    result, _, cached = lookup(cache, tmp_path, f">P1\n{SEQ_A}\n")
    assert cached == nucimport_output(("P1", VALUES_A))
    assert "1 sequences have different results in different store files" in result.stderr
