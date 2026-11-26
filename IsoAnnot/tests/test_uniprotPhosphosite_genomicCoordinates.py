# Synthetic gene on chromosome 1 (+) with three coding exons of four codons: 101-112, 201-212 and 301-312.
# The canonical protein P0ISO1 (MASKTTTTGSYW) uses the three; isoform P0ISO1-2 (MASKGSYW) skips the second.
CANONICAL, ISOFORM = "ENSP00000000101", "ENSP00000000102"
UNIPROT_FASTA = ">sp|P0ISO1|ISO_HUMAN Synthetic\nMASKTTTTGSYW\n>sp|P0ISO1-2|ISO_HUMAN Isoform 2 of Synthetic\nMASKGSYW\n"
PARSED = ("P0ISO1\t9606\tDR\tEnsembl\tENST00000000101.1\tENSP00000000101.1\tENSG00000000101.1. [P0ISO1-1]\n"
          "P0ISO1\t9606\tDR\tEnsembl\tENST00000000102.1\tENSP00000000102.1\tENSG00000000101.1. [P0ISO1-2]\n"
          "P0ISO1\t9606\tFT\tMOD_RES\t2\t3\tPhosphoserine\n"       # serine 3, in the first exon
          "P0ISO1\t9606\tFT\tMOD_RES\t4\t5\tPhosphothreonine\n"    # threonine 5, in the skipped exon
          "P0ISO1\t9606\tFT\tMOD_RES\t9\t10\tPhosphoserine\n"      # serine 10, serine 6 of the isoform
          "P0ISO1\t9606\tFT\tREGION\t2\t10\tDisordered\n")         # goes through the skipped exon
PHOSPHOSITE = "header\n" * 4 + "\t".join(["ISO", "ISO", "P0ISO1", "ISO", "Y11-p"]) + "\n"


def cds_lines(protein_id, exons):
    return "".join(f'1\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tgene_id "ENSG00000000101"; protein_id "{protein_id}";\n'
                   for start, end in exons)


def run_genomic_coordinates(run_script, tmp_path, gtf_text, ensembl_fasta):
    """Runs the script on the synthetic gene and returns its log and the rows of the domain file."""
    files = {"uniprot.fasta": UNIPROT_FASTA, "pep.fa": ensembl_fasta, "parsed.txt": PARSED,
             "phosphosite.txt": PHOSPHOSITE, "ensembl.gtf": gtf_text}
    for name, text in files.items():
        (tmp_path / name).write_text(text)
    result = run_script("uniprotPhosphosite_genomicCoordinates.py", "--uniprot_fasta", tmp_path / "uniprot.fasta",
                        "--ensembl_fasta", tmp_path / "pep.fa", "--uniprot_parsed", tmp_path / "parsed.txt",
                        "--phosphosite_files", tmp_path / "phosphosite.txt", "--ensembl_gtf", tmp_path / "ensembl.gtf",
                        "--biomart_host", "https://www.ensembl.org", "--output_protein", tmp_path / "protein.txt",
                        "--output_domain", tmp_path / "domain.txt")
    assert result.returncode == 0, result.stderr
    rows = [line.split("\t") for line in (tmp_path / "domain.txt").read_text().splitlines()]
    return result.stderr, rows


def test_alternative_isoform_gets_the_features_at_the_codons_of_the_canonical_protein(run_script, tmp_path):
    """
    UniProt numbers the features on the canonical sequence. The isoform used to get them at the same
    residue number: serine 10 was lost, and a feature could land on another residue of the same kind.
    """
    gtf = cds_lines(CANONICAL, [(101, 112), (201, 212), (301, 312)]) + cds_lines(ISOFORM, [(101, 112), (301, 312)])
    fasta = f">{CANONICAL}.1 pep\nMASKTTTTGSYW\n>{ISOFORM}.1 pep\nMASKGSYW\n"
    log, rows = run_genomic_coordinates(run_script, tmp_path, gtf, fasta)
    places = {protein: [(row[1], row[4], row[5], row[10]) for row in rows if row[13] == protein]
              for protein in (CANONICAL, ISOFORM)}
    assert places[CANONICAL] == [("Phosphoserine", "107", "109", "S"), ("Phosphothreonine", "201", "203", "T"),
                                 ("Phosphoserine", "304", "306", "S"), ("Disordered", "107", "306", "SKTTTTGS"),
                                 ("phosphotyrosine", "307", "309", "Y")]
    # The threonine and the region need the skipped exon, so the isoform does not get them
    assert places[ISOFORM] == [("Phosphoserine", "107", "109", "S"), ("Phosphoserine", "304", "306", "S"),
                               ("phosphotyrosine", "307", "309", "Y")]
    assert "match an alternative UniProt isoform" not in log


def test_alternative_isoform_in_another_reading_frame_gets_no_feature(run_script, tmp_path):
    """The isoform starts one base later in the last exon, so the codon of serine 10 is not a codon in it."""
    gtf = cds_lines(CANONICAL, [(101, 112), (201, 212), (301, 312)]) + cds_lines(ISOFORM, [(101, 112), (302, 313)])
    fasta = f">{CANONICAL}.1 pep\nMASKTTTTGSYW\n>{ISOFORM}.1 pep\nMASKGSYW\n"
    log, rows = run_genomic_coordinates(run_script, tmp_path, gtf, fasta)
    assert [(row[1], row[4], row[5]) for row in rows if row[13] == ISOFORM] == [("Phosphoserine", "107", "109")]


def test_isoform_without_a_canonical_protein_gets_no_feature_and_is_reported(run_script, tmp_path):
    """Without a protein of the canonical sequence, nothing says where the features are on the genome."""
    gtf = cds_lines(ISOFORM, [(101, 112), (301, 312)])
    fasta = f">{ISOFORM}.1 pep\nMASKGSYW\n"
    log, rows = run_genomic_coordinates(run_script, tmp_path, gtf, fasta)
    assert rows == []
    assert f"1 proteins match an alternative UniProt isoform" in log and ISOFORM in log


def test_stops_when_the_table_does_not_know_the_chromosomes_of_the_refseq_gtf(run_script, synthetic_dir, tmp_path):
    """RefSeq CDSs on sequences that chr2acc can't translate were kept under a name the genome does not have."""
    files = {"uniprot.fasta": UNIPROT_FASTA, "pep.fa": f">{CANONICAL}.1 pep\nMASKTTTTGSYW\n", "parsed.txt": PARSED,
             "phosphosite.txt": PHOSPHOSITE, "ensembl.gtf": cds_lines(CANONICAL, [(101, 112), (201, 212), (301, 312)]),
             "refseq.gtf": cds_lines("NP_900101.1", [(101, 112)]).replace("1\ttest", "NW_999999.1\ttest")}
    for name, text in files.items():
        (tmp_path / name).write_text(text)
    result = run_script("uniprotPhosphosite_genomicCoordinates.py", "--uniprot_fasta", tmp_path / "uniprot.fasta",
                        "--ensembl_fasta", tmp_path / "pep.fa", "--uniprot_parsed", tmp_path / "parsed.txt",
                        "--phosphosite_files", tmp_path / "phosphosite.txt", "--ensembl_gtf", tmp_path / "ensembl.gtf",
                        "--refseq_gtf", tmp_path / "refseq.gtf", "--chr_ref", synthetic_dir / "chr2acc",
                        "--biomart_host", "e115.ensembl.org", "--output_protein", tmp_path / "protein.txt",
                        "--output_domain", tmp_path / "domain.txt")
    assert result.returncode != 0
    assert "Critical Error in uniprotPhosphosite_genomicCoordinates.py: 1/1 chromosome lookups failed" in result.stderr
