import os
import sys

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)


def read_fasta(path):
    seqs = {}
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line)
    return {name: "".join(lines) for name, lines in seqs.items()}


def test_refseq_transcripts_are_spliced_from_the_genome(run_script, synthetic_dir, tmp_path):
    """
    The synthetic rna.fna is the spliced exons, so the output must be identical to it, on both strands.
    The GTF has RefSeq chromosome names (NC_000001.11) and the genome Ensembl ones (1).
    """
    output = tmp_path / "refseq_transcripts.fa"
    result = run_script("transcripts_from_genome.py", "--gtf", synthetic_dir / "refseq.gtf",
                        "--genome", synthetic_dir / "genome.fa", "--chr_ref", synthetic_dir / "chr2acc",
                        "--output", output)
    assert result.returncode == 0, result.stderr
    assert read_fasta(output) == read_fasta(synthetic_dir / "refseq_rna.fna")
    assert max(len(line) for line in output.read_text().splitlines()) == 60
    assert not os.path.exists(f"{output}.tmp")


def test_transcripts_off_the_genome_are_left_out(run_script, synthetic_dir, tmp_path):
    """A transcript on a chromosome the genome does not have is skipped with one warning, not written empty."""
    gtf = tmp_path / "refseq.gtf"
    gtf.write_text((synthetic_dir / "refseq.gtf").read_text()
                   + 'NC_000002.12\tBestRefSeq\texon\t1\t100\t.\t+\t.\tgene_id "OTHER"; transcript_id "NM_999999.1";\n')
    output = tmp_path / "refseq_transcripts.fa"
    result = run_script("transcripts_from_genome.py", "--gtf", gtf, "--genome", synthetic_dir / "genome.fa",
                        "--chr_ref", synthetic_dir / "chr2acc", "--output", output)
    assert result.returncode == 0, result.stderr
    assert sorted(read_fasta(output)) == ["NM_900001.1", "NM_900002.1", "NM_900003.2", "NM_900004.1",
                                          "NM_900006.3", "NR_900005.1"]
    assert "1 transcripts are on chromosomes that are not in the genome (NC_000002.12)" in result.stderr
