import os
import shutil
import subprocess
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

ISOANNOT_DIR = os.path.dirname(main_path)
# ISOANNOT_SNAKEMAKE, or snakemake on PATH
SNAKEMAKE = os.environ.get("ISOANNOT_SNAKEMAKE") or shutil.which("snakemake")

# First two lines of every NucImport output, as NucImportMay2012.jar writes them
NUCIMPORT_HEADER = ("Protein ID\tImport Probability\t cNLS probability\tPredicted Most probable cNLS\t"
                    "Predicted most probable cNLS class number\tPosition")
NUCIMPORT_STARS = "*" * 93


@pytest.fixture
def nls_dir(tmp_path):
    """nls folder of the ensembl hsapiens run whose path_output is tmp_path."""
    return tmp_path / "data" / "Hsapiens" / "output" / "ensembl" / "nls"


@pytest.fixture
def run_rule(tmp_path):
    """
    Runs the Snakefile rules for the given targets, with the ensembl hsapiens config and path_output=tmp_path.
    Only the given rule is allowed, so its input files are used as they are and nothing upstream runs.
    Scripts run with the test interpreter (first python3 on PATH). Returns the CompletedProcess.
    """
    if not SNAKEMAKE or not os.path.isfile(SNAKEMAKE):
        pytest.skip("Snakemake not found (set ISOANNOT_SNAKEMAKE or put snakemake on PATH)")

    def _run(rule, targets, *config, dry_run=False):
        config_dir = os.path.join(ISOANNOT_DIR, "config", "ensembl", "hsapiens")
        cmd = [SNAKEMAKE, "-p",
               "--snakefile", os.path.join(config_dir, "Snakefile.smk"),
               "--configfile", os.path.join(config_dir, "config.yaml"),
               "--config", "db=ensembl", f"path_output={tmp_path}", *config,
               "--directory", ISOANNOT_DIR, "--cores", "2", "--allowed-rules", rule,
               "--nolock", "--drop-metadata", *[str(t) for t in targets]]
        if dry_run:
            cmd.insert(1, "-n")
        env = dict(os.environ, PATH=os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"])
        return subprocess.run(cmd, capture_output=True, text=True, env=env)
    return _run


def nucimport_row(prot_id):
    return f"{prot_id}\t0.9949979710297792\t0.806466703423478\tAALLLLKRKRRSPPKKKRKV\tClass1\t15\t"


def nucimport_output(*prot_ids):
    """NucImport output: a header line, a line of stars, then one line per protein."""
    return "\n".join([NUCIMPORT_HEADER, NUCIMPORT_STARS] + [nucimport_row(p) for p in prot_ids]) + "\n"


def test_merge_nls_chunks_keeps_first_protein(run_rule, nls_dir):
    """The merge skipped 3 lines of every chunk after the first, so it lost their first protein."""
    chunk_dir = nls_dir / "tmp"
    chunk_dir.mkdir(parents=True)
    (nls_dir / "cached_results.txt").write_text(nucimport_output())
    (chunk_dir / "output_chunk_0.txt").write_text(nucimport_output("P1", "P4"))
    (chunk_dir / "output_chunk_1.txt").write_text(nucimport_output("P2", "P5"))
    (chunk_dir / "output_chunk_2.txt").write_text(nucimport_output("P3"))
    merged = nls_dir / "output_merged.txt"
    result = run_rule("merge_nls_chunks", [merged], "nls_chunks=3")
    assert result.returncode == 0, result.stdout + result.stderr
    assert merged.read_text() == nucimport_output("P1", "P4", "P2", "P5", "P3")


def test_merge_nls_chunks_with_cached(run_rule, nls_dir):
    """The cached results come first, then the new ones; with every chunk empty the header still comes from the cached file."""
    chunk_dir = nls_dir / "tmp"
    chunk_dir.mkdir(parents=True)
    (nls_dir / "cached_results.txt").write_text(nucimport_output("P1", "P2"))
    (chunk_dir / "output_chunk_0.txt").write_text(nucimport_output("P3"))
    (chunk_dir / "output_chunk_1.txt").write_text("")
    merged = nls_dir / "output_merged.txt"
    result = run_rule("merge_nls_chunks", [merged], "nls_chunks=2")
    assert result.returncode == 0, result.stdout + result.stderr
    assert merged.read_text() == nucimport_output("P1", "P2", "P3")


def test_split_proteins_fewer_proteins_than_chunks(run_rule, nls_dir):
    """With fewer proteins than nls_chunks, awk did not create the extra chunks and Snakemake stopped with MissingOutputException."""
    nls_dir.mkdir(parents=True)
    (nls_dir / "to_predict.fa").write_text(">P1\nMSPKKKRKVEDPAAAALLLL\n>P2\nMAAAAGGGGLLLLVVVVSSSS\n")
    chunks = [nls_dir / "chunks_temp" / f"chunk_{n}.fa" for n in range(4)]
    result = run_rule("split_proteins", chunks, "nls_chunks=4")
    assert result.returncode == 0, result.stdout + result.stderr
    records = [[line for line in chunk.read_text().splitlines() if line] for chunk in chunks]
    assert records == [[">P1", "MSPKKKRKVEDPAAAALLLL"], [">P2", "MAAAAGGGGLLLLVVVVSSSS"], [], []]


@pytest.mark.parametrize("empty_output", ["", nucimport_output()])
def test_merge_nls_chunks_empty_chunks(run_rule, nls_dir, empty_output):
    """
    Empty chunks give an empty output (run_nucimport skips Java), or only the two header lines when NucImport
    ran on them. The pattern-based merge keeps one header in both cases.
    """
    chunk_dir = nls_dir / "tmp"
    chunk_dir.mkdir(parents=True)
    (nls_dir / "cached_results.txt").write_text(nucimport_output())
    (chunk_dir / "output_chunk_0.txt").write_text(empty_output)
    (chunk_dir / "output_chunk_1.txt").write_text(nucimport_output("P1"))
    (chunk_dir / "output_chunk_2.txt").write_text(empty_output)
    merged = nls_dir / "output_merged.txt"
    result = run_rule("merge_nls_chunks", [merged], "nls_chunks=3")
    assert result.returncode == 0, result.stdout + result.stderr
    assert merged.read_text() == nucimport_output("P1")


def folder_state(path):
    return sorted((entry.name, entry.stat().st_mtime_ns) for entry in os.scandir(path))


@pytest.fixture
def leftover_aout():
    """
    software/NucImport/aout.txt as left by the old run_nucimport, which ran the jar in that folder.
    It is untracked, so installs that ran NLS before keep it after a pull. A real leftover is used as it is;
    otherwise one is written and removed after the test.
    """
    aout = os.path.join(ISOANNOT_DIR, "software", "NucImport", "aout.txt")
    if os.path.exists(aout):
        yield aout
        return
    with open(aout, "w") as h:
        h.write("OLD RUN LEFTOVER\n")
    yield aout
    os.remove(aout)


def test_run_nucimport_own_folder(run_rule, nls_dir, leftover_aout):
    """Parallel NucImport runs all wrote aout.txt into software/NucImport; each chunk now runs in its own folder, also when an old aout.txt is there. Empty chunks give an empty output."""
    if shutil.which("java") is None:
        pytest.skip("java not found on PATH")
    chunk_dir = nls_dir / "chunks_temp"
    chunk_dir.mkdir(parents=True)
    (chunk_dir / "chunk_0.fa").write_text(">P1\nMSPKKKRKVEDPAAAALLLLKRKRRSPPKKKRKVAAAGGGLLL\n"
                                          ">P3\nMKRPRGRPKKKRKVLLAAAPEDKRRKAAAGGGLLLAAAGGGLLLDDE\n")
    (chunk_dir / "chunk_1.fa").write_text("")
    nucimport_dir = os.path.join(ISOANNOT_DIR, "software", "NucImport")
    before = folder_state(nucimport_dir)
    with open(leftover_aout) as h:
        aout_before = h.read()
    outputs = [nls_dir / "tmp" / f"output_chunk_{n}.txt" for n in range(2)]
    result = run_rule("run_nucimport", outputs, "nls_chunks=2")
    assert result.returncode == 0, result.stdout + result.stderr
    assert folder_state(nucimport_dir) == before
    with open(leftover_aout) as h:
        assert h.read() == aout_before
    assert sorted(os.listdir(nls_dir / "tmp")) == ["output_chunk_0.txt", "output_chunk_1.txt"]
    lines = outputs[0].read_text().splitlines()
    assert lines[:2] == [NUCIMPORT_HEADER, NUCIMPORT_STARS]
    rows = [line.split("\t") for line in lines[2:]]
    assert [[r[0], r[3], r[4], r[5]] for r in rows] == [["P1", "AALLLLKRKRRSPPKKKRKV", "Class1", "15"],
                                                         ["P3", "PRGRPKKKRKVLLAAAPEDK", "Class1", "4"]]
    assert [float(r[1]) for r in rows] == pytest.approx([0.99499797, 0.99333969])
    assert outputs[1].read_text() == ""


def test_run_nucimport_empty_chunk(run_rule, tmp_path, nls_dir):
    """An empty chunk gets an empty output without starting Java, which took about 6 s per chunk."""
    chunk_dir = nls_dir / "chunks_temp"
    chunk_dir.mkdir(parents=True)
    (chunk_dir / "chunk_0.fa").write_text("")
    output = nls_dir / "tmp" / "output_chunk_0.txt"
    result = run_rule("run_nucimport", [output], "nls_chunks=1")
    assert result.returncode == 0, result.stdout + result.stderr
    assert output.read_text() == ""
    log = tmp_path / "logs" / "Hsapiens" / "ensembl" / "run_nucimport" / "run_nucimport_chunk_0.log"
    assert log.read_text() == "Empty chunk: NucImport not run\n"
    assert not (nls_dir / "tmp" / "nucimport_chunk_0").exists()


@pytest.mark.parametrize("config, xmx", [((), "-Xmx2g"), (("nls_java_memory=6g",), "-Xmx6g")])
def test_run_nucimport_java_memory(run_rule, nls_dir, config, xmx):
    """Each NucImport JVM gets a maximum heap, 2g by default as the NucImport README asks."""
    chunk_dir = nls_dir / "chunks_temp"
    chunk_dir.mkdir(parents=True)
    (chunk_dir / "chunk_0.fa").write_text(">P1\nMSPKKKRKVEDPAAAALLLLKRKRRSPPKKKRKVAAAGGGLLL\n")
    result = run_rule("run_nucimport", [nls_dir / "tmp" / "output_chunk_0.txt"], "nls_chunks=1", *config,
                      dry_run=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"java {xmx} -jar " in result.stdout


def test_expand_nls_writes_log(run_rule, tmp_path, nls_dir):
    """expand_nls declared a log but never wrote to it."""
    nls_dir.mkdir(parents=True)
    (nls_dir / "nls_parsed.tsv").write_text("ID\tImpProb\tcNLSProb\tSeq\tClass\tPos\n"
                                            "P1\t0.995\t0.806\tAALLLLKRKRRSPPKKKRKV\tClass1\t15\n")
    (nls_dir / "protein_mapping.tsv").write_text("UniqueID\tOriginalIDs\nP1\tP1,P2\nP3\tP3\n")
    (nls_dir / "nucimport_missing.txt").write_text("")
    result = run_rule("expand_nls", [nls_dir / "nls_final_expanded.tsv"])
    assert result.returncode == 0, result.stdout + result.stderr
    log = tmp_path / "logs" / "Hsapiens" / "ensembl" / "expand_nls.log"
    assert log.read_text().splitlines()[-1] == "INFO: Expansion complete. Total annotations generated: 2"


@pytest.mark.parametrize("merged_text", ["", nucimport_output()])
def test_parse_nls_no_results(run_rule, tmp_path, nls_dir, merged_text):
    """No protein reached NucImport (all shorter than 30 aa): parse_nls stopped the run instead of giving an empty layer."""
    nls_dir.mkdir(parents=True)
    (nls_dir / "output_merged.txt").write_text(merged_text)
    parsed = nls_dir / "nls_parsed.tsv"
    result = run_rule("parse_nls", [parsed])
    assert result.returncode == 0, result.stdout + result.stderr
    assert parsed.read_text() == "ID\tImpProb\tcNLSProb\tSeq\tClass\tPos\n"
    log = (tmp_path / "logs" / "Hsapiens" / "ensembl" / "parse_nls.log").read_text()
    assert "WARNING: NucImport returned no results, so the NLS layer will be empty." in log


def test_expand_nls_no_results(run_rule, nls_dir):
    """The empty table of parse_nls gives an empty expanded table."""
    nls_dir.mkdir(parents=True)
    (nls_dir / "nls_parsed.tsv").write_text("ID\tImpProb\tcNLSProb\tSeq\tClass\tPos\n")
    (nls_dir / "protein_mapping.tsv").write_text("UniqueID\tOriginalIDs\nP1\tP1,P2\n")
    (nls_dir / "nucimport_missing.txt").write_text("")
    result = run_rule("expand_nls", [nls_dir / "nls_final_expanded.tsv"])
    assert result.returncode == 0, result.stdout + result.stderr
    assert (nls_dir / "nls_final_expanded.tsv").read_text() == "ID\tImpProb\tcNLSProb\tSeq\tClass\tPos\n"


def test_nls_deduplicate_no_proteins(run_rule, nls_dir):
    """No protein of 30 aa or longer: the redundancy percentage divided by zero."""
    nls_dir.mkdir(parents=True)
    (nls_dir / "nls_filtered_proteins.fa").write_text("")
    result = run_rule("nls_deduplicate", [nls_dir / "unique_proteins.fa"])
    assert result.returncode == 0, result.stdout + result.stderr
    assert (nls_dir / "unique_proteins.fa").read_text() == ""
    assert (nls_dir / "protein_mapping.tsv").read_text() == "UniqueID\tOriginalIDs\n"
