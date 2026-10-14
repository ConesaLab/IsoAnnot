# Synthetic test dataset

This folder holds a small synthetic dataset for the IsoAnnot tests. It has one 20 kb chromosome,
four genes and six transcripts, written the way Ensembl, RefSeq, NCBI and UniProt write their files.
It is small enough to run any IsoAnnot script in a fraction of a second, and every transcript is there to
exercise a specific case.

All files are written by `make_fixture.py` with a fixed seed, so they are the same on every run. If you
change the transcript models, or one of the scripts that make the derived files (`referenceSQANTI.py`,
`parse_protein_databases.py`), run it again and commit the result:

```bash
python3 IsoAnnot/tests/fixtures/synthetic/make_fixture.py
```

`test_synthetic_fixture.py::test_fixture_is_up_to_date` fails when the committed files no longer match
what the generator writes.

In the tests, use the `synthetic_dir` fixture to get this folder and `run_script` to run a script from
`IsoAnnot/scripts/` (both in `tests/conftest.py`). Treat the folder as read-only and write outputs to
`tmp_path`.

## Transcripts

All coordinates are 1-based and inclusive. The CDS does not include the stop codon, as in the Ensembl and
RefSeq GTFs.

| Ensembl transcript | RefSeq | Gene | Strand | Exons | CDS | Protein (Ensembl / RefSeq) | UniProt | Case |
|---|---|---|---|---|---|---|---|---|
| `ENST00000000001` (`SYNA-201`) | `NM_900001.1` | `SYNA` | + | 3 | 1101-3200 | `ENSP00000000001.2` / `NP_900001.1` | `P0SYA1` (SwissProt canonical) | Multi-exon coding, + strand |
| `ENST00000000002` (`SYNA-202`) | `NM_900002.1` | `SYNA` | + | 4 | 1101-2511 | `ENSP00000000002.1` / `NP_900002.1` | `A0A0S1YNA6` (TrEMBL) | NMD: the extra exon 2501-2600 shifts the frame and the stop codon is 86 nt before the last junction |
| `ENST00000000003` (`SYNB-201`) | `NM_900003.2` | `SYNB` | - | 4 | 6100-9650 | `ENSP00000000003.2` / `NP_900003.1` | `P0SYB1-2` (SwissProt isoform) | Multi-exon coding, - strand |
| `ENST00000000004` (`SYNB-202`) | `NM_900004.1` | `SYNB` | - | 5 | 6100-9650 | `ENSP00000000004.1` / `NP_900004.1` | `P0SYB1-2` (SwissProt isoform) | Same CDS as `SYNB-201` with a different 5'UTR: two protein IDs with an identical sequence |
| `ENST00000000005` (`SYNC-201`) | `NR_900005.1` | `SYNC` | + | 2 | none | none | none | Non-coding (`lncRNA`) |
| `ENST00000000006` (`SYND-201`) | `NM_900006.3` | `SYND` | + | 1 | 14601-15401 | `ENSP00000000006.3` / `NP_900006.2` | none | Single-exon coding, no UniProt match |

Other details:
* Every intron has canonical `GT..AG` splice sites (13 introns in total).
* Exons shared by several transcripts have the same Ensembl `exon_id`, as in a real Ensembl GTF.
* The canonical SwissProt entry of SYNB, `P0SYB1`, is the `P0SYB1-2` sequence with 10 extra residues, so it
  matches no transcript.
* IDs, versions and accessions are made up. The `NM_9xxxxx`, `NP_9xxxxx`, `P0SYxx` and `A0A0S1YNA6`
  accessions do not exist in RefSeq or UniProt.

## Input files

These files stand for what the `get_*` and `prepare_*` rules leave in `data/<prefix>/config/`. They are
plain text: most real downloads are gzipped, and the pipeline decompresses them before use.

* **`genome.fa`**: Ensembl genome FASTA with one chromosome, `1`, of 20,000 bp (60 bases per line).
  Header: `>1 dna:chromosome chromosome:SYNTH1:1:1:20000:1 REF`. The background is random; the exons, ORFs
  and splice sites are written on top of it.
* **`ensembl.gtf`**: Ensembl GTF with `#!` header lines and `gene`, `transcript`, `exon`, `CDS`,
  `start_codon`, `stop_codon`, `five_prime_utr` and `three_prime_utr` lines. IDs are unversioned, with the
  version in `gene_version`, `transcript_version` and `protein_version`, as Ensembl writes them.
* **`ensembl_cdna.fa`**: Ensembl `cdna.all` style FASTA, one sequence per transcript, with versioned
  headers (`>ENST00000000001.4 cdna chromosome:SYNTH1:1:1001:3400:1 gene:ENSG00000000001.3 ...`).
  The real `cdna.all` file leaves lncRNAs out (they are in `ncrna`); here `SYNC-201` is included so the same
  file can be used as a `mytranscripts` input (`--config fasta_cdna=...`).
* **`ensembl_pep.fa`**: Ensembl `pep.all` style FASTA with versioned headers (`>ENSP00000000001.2 pep ...`).
  Each sequence is the translation of the CDS, without the stop `*`.
* **`refseq.gtf`**: the same models in NCBI RefSeq GTF format on `NC_000001.11`, as `prepare_refseq_gtf`
  leaves it (only `NC_` lines, no header). `gene_id` is the gene symbol, the NCBI gene ID is in
  `db_xref "GeneID:..."`, and transcript and protein IDs are versioned (`NM_900001.1`, `NP_900001.1`).
* **`refseq_rna.fna`**: RefSeq `rna.fna` style FASTA (`>NM_900001.1 Homo sapiens synthetic gene A (SYNA), ...`).
* **`refseq_protein.faa`**: RefSeq `protein.faa` style FASTA (`>NP_900001.1 synthetic protein A isoform a [Homo sapiens]`).
* **`chr2acc`**: NCBI chromosome to accession table, tab-separated, with its `#Chromosome` /
  `Accession.version` header line. It maps `1` to `NC_000001.11`.
* **`uniprot.fasta`** and **`uniprot_additional.fasta`**: UniProt reference proteome FASTA files. The first
  one has the SwissProt canonical entries (`sp|...`) and one TrEMBL entry (`tr|...`), the second one the
  SwissProt isoform `sp|P0SYB1-2|...`, as in the `UP..._9606.fasta` and `UP..._9606_additional.fasta` pair.

## Derived files

These files are made by running the pipeline scripts on the input files, with the same arguments as the
Snakefile rules. They are the exact output of the current scripts, including their known problems (for
example, `ORF_length` is written as a float).

* **`ensembl_sqanti_classification.txt`**, **`ensembl_sqanti_junctions.txt`** and **`ensembl_sqanti_NMD.txt`**:
  output of `referenceSQANTI.py --database ensembl` (rule `run_refsqanti`). `refProt` is the unversioned
  Ensembl protein ID, as in a real run. `ENST00000000002` is the only `NMD` transcript.
* **`refseq_sqanti_classification.txt`**, **`refseq_sqanti_junctions.txt`** and **`refseq_sqanti_NMD.txt`**:
  output of `referenceSQANTI.py --database refseq`. The `NC_000001.11` chromosome is written as `1`, the
  name in the Ensembl genome. `refProt` is the versioned RefSeq protein ID and `structural_category`
  is `NA`. `NM_900002.1` is the only `NMD` transcript.
* **`parsed_databases.json`**: output of `parse_protein_databases.py` (rule `parse_protein_databases`). It maps
  each protein sequence to its IDs by source (`swissprot`, `trembl`, `refseq`, `ensembl`), without versions.

## Not included

The dataset does not have SQANTI3 outputs for `mytranscripts` (corrected GTF and ORF `.faa`), UniProt `.dat`
files, InterProScan, UTRscan, RepeatMasker, NucImport or miRWalk outputs. Tests that need them should build
small ones inline with `tmp_path`, or add them to the generator if several tests share them.
