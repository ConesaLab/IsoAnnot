import os
import subprocess
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from repeatmasker_ids import shorten_fasta_ids, restore_out_ids

ISOANNOT_DIR = os.path.dirname(main_path)
SNAKEMAKE = os.environ.get("ISOANNOT_SNAKEMAKE", "/media/pabloati/HDD1/conda/envs/snakemake/bin/snakemake")

LONG_ID = "PB.1.1|chr1:1000-5000(+)|transcript/1234|full_length_coverage=10|length=4000"
SHORT_ID = "ENST00000000001.4"

# Stand-in for RepeatMasker: it stops on IDs longer than 50 characters, as the real one does, and writes
# one .out hit per sequence, with the 3 header lines and the column layout of a real .out file.
FAKE_REPEATMASKER = f"""#!{sys.executable}
import os, sys
fasta = sys.argv[1]
out_dir = sys.argv[sys.argv.index("-dir") + 1]
ids = [line[1:].split()[0] for line in open(fasta) if line.startswith(">")]
if any(len(i) > 50 for i in ids):
    sys.exit("FastaDB::_cleanIDs: Fasta file contains a sequence identifier which is too long ( max id length = 50 )")
with open(os.path.join(out_dir, os.path.basename(fasta) + ".out"), "w") as out:
    out.write("   SW   perc perc perc  query      position in query    matching  repeat        position in repeat\\n")
    out.write("score   div. del. ins.  sequence   begin  end  (left)   repeat    class/family  begin  end (left)  ID\\n")
    out.write("\\n")
    for n, i in enumerate(ids, 1):
        out.write(f"  239   29.4  1.9  1.0  {{i}}  11  160  (840) +  AluY  SINE/Alu  1  148  (163)  {{n}}\\n")
"""

OUT_HEADER = ("   SW   perc perc perc  query      position in query    matching  repeat        position in repeat\n"
              "score   div. del. ins.  sequence   begin  end  (left)   repeat    class/family  begin  end (left)  ID\n"
              "\n")


def fake_out_line(seq_id, n):
    return f"  239   29.4  1.9  1.0  {seq_id}  11  160  (840) +  AluY  SINE/Alu  1  148  (163)  {n}\n"


def test_round_trip(tmp_path):
    """IDs over 50 characters are hashed for RepeatMasker and restored in its .out file."""
    id_50 = "X" * 50
    fasta = tmp_path / "cdna.fa"
    fasta.write_text(f">{SHORT_ID} cdna\nACGT\n>{LONG_ID} extra words\nGGCC\nTTAA\n>{id_50}\nCCCC\n")
    hashed, mapping = tmp_path / "tmp_hashed.fasta", tmp_path / "id_mapping.tsv"

    assert shorten_fasta_ids(str(fasta), str(hashed), str(mapping)) == 1
    short_id = "dced6ef2571fd5ad37a77374b1ba204a"  # MD5 of LONG_ID
    assert mapping.read_text() == f"{short_id}\t{LONG_ID}\n"
    assert hashed.read_text() == f">{SHORT_ID} cdna\nACGT\n>{short_id}\nGGCC\nTTAA\n>{id_50}\nCCCC\n"

    rm_out = tmp_path / "tmp_hashed.fasta.out"
    rm_out.write_text(OUT_HEADER + fake_out_line(SHORT_ID, 1) + fake_out_line(short_id, 2) + fake_out_line(id_50, 3))
    restored = tmp_path / "cdna.fa.out"
    assert restore_out_ids(str(rm_out), str(mapping), str(restored)) == 1
    assert restored.read_text() == OUT_HEADER + fake_out_line(SHORT_ID, 1) + fake_out_line(LONG_ID, 2) + fake_out_line(id_50, 3)
    assert not os.path.exists(f"{restored}.tmp")


def test_restored_out_parsed_by_layer(tmp_path):
    """layer_repeatmasker.py reads the restored .out with the original ID in the transcript column."""
    from layer_repeatmasker import get_repeatmasker_info
    mapping = tmp_path / "id_mapping.tsv"
    mapping.write_text(f"0123456789abcdef0123456789abcdef\t{LONG_ID}\n")
    rm_out = tmp_path / "tmp_hashed.fasta.out"
    rm_out.write_text(OUT_HEADER + fake_out_line("0123456789abcdef0123456789abcdef", 1))
    restored = tmp_path / "cdna.fa.out"
    restore_out_ids(str(rm_out), str(mapping), str(restored))

    assert list(get_repeatmasker_info(str(restored))) == [
        {"identifier": "SINE/Alu", "transcript": LONG_ID, "start": 11, "end": 160, "desc": "AluY"}]


def test_short_ids_and_no_repeats(tmp_path, synthetic_dir):
    """With no long IDs the FASTA is copied unchanged, and a .out without hits is copied as is."""
    hashed, mapping = tmp_path / "tmp_hashed.fasta", tmp_path / "id_mapping.tsv"
    assert shorten_fasta_ids(str(synthetic_dir / "ensembl_cdna.fa"), str(hashed), str(mapping)) == 0
    assert hashed.read_text() == (synthetic_dir / "ensembl_cdna.fa").read_text()
    assert mapping.read_text() == ""

    rm_out = tmp_path / "tmp_hashed.fasta.out"
    rm_out.write_text("There were no repetitive sequences detected in tmp_hashed.fasta\n")
    restored = tmp_path / "cdna.fa.out"
    assert restore_out_ids(str(rm_out), str(mapping), str(restored)) == 0
    assert restored.read_text() == rm_out.read_text()


def test_run_repeatmasker_rule(tmp_path):
    """run_repeatmasker runs end to end and the .out file keeps the original long IDs."""
    if not os.path.isfile(SNAKEMAKE):
        pytest.skip(f"Snakemake not found at {SNAKEMAKE} (set ISOANNOT_SNAKEMAKE to use another one)")

    out_dir = tmp_path / "out"
    cdna = out_dir / "data" / "Hsapiens" / "config" / "ensembl" / "Homo_sapiens.GRCh38.cdna.all.fa"
    cdna.parent.mkdir(parents=True)
    cdna.write_text(f">{SHORT_ID} cdna\nACGTACGT\n>{LONG_ID}\nGGGGCCCC\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "RepeatMasker").write_text(FAKE_REPEATMASKER)
    (bin_dir / "RepeatMasker").chmod(0o755)
    # The rule links libnsl into $CONDA_PREFIX/lib; point it at a scratch folder instead of a real env.
    (tmp_path / "conda" / "lib").mkdir(parents=True)
    env = dict(os.environ, CONDA_PREFIX=str(tmp_path / "conda"),
               PATH=f"{bin_dir}:{os.path.dirname(sys.executable)}:{os.environ['PATH']}")
    # Conda activation does not set LD_LIBRARY_PATH, and the rule must work without it
    env.pop("LD_LIBRARY_PATH", None)

    config_dir = os.path.join(ISOANNOT_DIR, "config", "ensembl", "hsapiens")
    target = out_dir / "data" / "Hsapiens" / "output" / "ensembl" / "repeat_masker" / f"{cdna.name}.out"
    cmd = [SNAKEMAKE, "--snakefile", os.path.join(config_dir, "Snakefile.smk"),
           "--configfile", os.path.join(config_dir, "config.yaml"),
           "--config", "db=ensembl", f"path_output={out_dir}",
           "--directory", ISOANNOT_DIR, "--cores", "1", "--nolock", str(target)]
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    log = out_dir / "logs" / "Hsapiens" / "ensembl" / "run_repeatmasker.log"
    assert result.returncode == 0, result.stderr + (log.read_text() if log.exists() else "")

    rows = [line.split() for line in target.read_text().splitlines()[3:]]
    assert [row[4] for row in rows] == [SHORT_ID, LONG_ID]
