# A transcript of 30 bases in two exons. NCBI places it on chromosome 1 and, as NM_900010.2_1, on an
# alternate locus, where the intron is 100 bases longer.
SEQUENCE = "ACGTTGCAAGGCTTAACCGGATCGATTACG"
GTF = "".join(f'{chrom}\tBestRefSeq\texon\t{start}\t{end}\t.\t+\t.\tgene_id "SYNM"; transcript_id "{transcript}";\n'
              for chrom, start, end, transcript in [("NC_000001.11", 1001, 1010, "NM_900010.2"), ("NC_000001.11", 1101, 1120, "NM_900010.2"),
                                                    ("NT_900001.1", 501, 510, "NM_900010.2_1"), ("NT_900001.1", 701, 720, "NM_900010.2_1")])
# miRWalk leaves out the first base of the site: positions 6-15 of the transcript, duplex target = bases 7-15
MIRWALK = ("miRNA\tmRNA\tGenesymbol\tbinding_site\tduplex\tscore\tsite_type\n"
           f"hsa-miR-1\tNM_900010\tSYNM\t6,15\tTGGAATGT#{SEQUENCE[6:15]}#((((((((#))))))))\t1\t8mer\n")


def test_site_is_placed_on_each_placement_of_the_transcript(run_script, synthetic_dir, tmp_path):
    """
    The exons of the two placements were merged into one model of 60 bases, so the site was placed with
    the wrong exons. Each placement now gets the site with its own exons.
    """
    for name, text in {"rna.fna": f">NM_900010.2 SYNM\n{SEQUENCE}\n", "genomic.gtf": GTF, "mirwalk.txt": MIRWALK}.items():
        (tmp_path / name).write_text(text)
    output = tmp_path / "interactions_gc.txt"
    result = run_script("mirna_bs_genomic_coord.py", "--mirwalk_file", tmp_path / "mirwalk.txt", "--refseq_fasta", tmp_path / "rna.fna",
                        "--refseq_gtf", tmp_path / "genomic.gtf", "--mirna_output", output, "--chr_ref", synthetic_dir / "chr2acc")
    assert result.returncode == 0, result.stderr
    rows = [line.split("\t")[2:9] for line in output.read_text().splitlines()[1:]]
    # Bases 6-10 are the end of the first exon and 11-15 the start of the second
    assert rows == [["NM_900010", "1", "+", "1006", "1105", "1006,1101", "1010,1105"],
                    ["NM_900010", "NT_900001.1", "+", "506", "705", "506,701", "510,705"]]
