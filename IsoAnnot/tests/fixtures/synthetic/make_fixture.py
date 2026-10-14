#!/usr/bin/env python3
"""
Writes the synthetic IsoAnnot test dataset: a 20 kb genome with one chromosome, four genes
and six transcripts in Ensembl and RefSeq form, their cDNAs and proteins, a few UniProt entries and
the NCBI chr2acc table. Then it runs referenceSQANTI.py and parse_protein_databases.py on them, with
the same arguments as the Snakefile rules, so the derived files have the exact format of a real run.

The output is the same on every run (fixed seed). Run it again after changing the transcript models
or the scripts that make the derived files, and commit the result:

    python3 IsoAnnot/tests/fixtures/synthetic/make_fixture.py
"""

import argparse, os, random, subprocess, sys

FIXTURE_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.abspath(os.path.join(FIXTURE_DIR, "..", "..", "..", "scripts"))

SEED = 2026
GENOME_LENGTH = 20000
CHROM = "1"
REFSEQ_CHROM = "NC_000001.11"
ASSEMBLY = "SYNTH1"

BASES = "TCAG"
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: AMINO_ACIDS[16 * i + 4 * j + k]
               for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
STOP_CODONS = sorted(c for c, aa in CODON_TABLE.items() if aa == "*")
SENSE_CODONS = sorted(c for c, aa in CODON_TABLE.items() if aa != "*")
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

# Transcript models. Coordinates are 1-based and inclusive, exons in ascending order.
# cds = (lowest, highest) genomic position of the CDS without the stop codon, as Ensembl and
# RefSeq write it. The stop codon is the next 3 bases downstream. See README.md for the cases.
GENES = [
    {
        "gene_id": "ENSG00000000001", "gene_version": 3, "name": "SYNA", "strand": "+",
        "biotype": "protein_coding", "ncbi_gene_id": 900001, "description": "synthetic gene A",
        "transcripts": [
            {"id": "ENST00000000001", "version": 4, "name": "SYNA-201", "biotype": "protein_coding",
             "tags": ["basic", "Ensembl_canonical"],
             "exons": [(1001, 1200), (2001, 2150), (3001, 3400)], "cds": (1101, 3200),
             "protein_id": "ENSP00000000001", "protein_version": 2,
             "refseq_id": "NM_900001.1", "refseq_protein": "NP_900001.1"},
            # Extra exon 2501-2600 shifts the frame: premature stop more than 55 nt before the
            # last exon-exon junction, so this transcript is NMD.
            {"id": "ENST00000000002", "version": 1, "name": "SYNA-202", "biotype": "nonsense_mediated_decay",
             "tags": [],
             "exons": [(1001, 1200), (2001, 2150), (2501, 2600), (3001, 3400)], "cds": (1101, 2511),
             "protein_id": "ENSP00000000002", "protein_version": 1,
             "refseq_id": "NM_900002.1", "refseq_protein": "NP_900002.1"},
        ],
    },
    {
        "gene_id": "ENSG00000000002", "gene_version": 5, "name": "SYNB", "strand": "-",
        "biotype": "protein_coding", "ncbi_gene_id": 900002, "description": "synthetic gene B",
        "transcripts": [
            {"id": "ENST00000000003", "version": 2, "name": "SYNB-201", "biotype": "protein_coding",
             "tags": ["basic", "Ensembl_canonical"],
             "exons": [(6001, 6300), (7001, 7120), (8001, 8090), (9501, 9700)], "cds": (6100, 9650),
             "protein_id": "ENSP00000000003", "protein_version": 2,
             "refseq_id": "NM_900003.2", "refseq_protein": "NP_900003.1"},
            # Same CDS as SYNB-201 with a different 5'UTR: two protein IDs, one identical sequence.
            {"id": "ENST00000000004", "version": 1, "name": "SYNB-202", "biotype": "protein_coding",
             "tags": ["basic"],
             "exons": [(6001, 6300), (7001, 7120), (8001, 8090), (9501, 9660), (10501, 10600)], "cds": (6100, 9650),
             "protein_id": "ENSP00000000004", "protein_version": 1,
             "refseq_id": "NM_900004.1", "refseq_protein": "NP_900004.1"},
        ],
    },
    {
        "gene_id": "ENSG00000000003", "gene_version": 1, "name": "SYNC", "strand": "+",
        "biotype": "lncRNA", "ncbi_gene_id": 900003, "description": "synthetic gene C",
        "transcripts": [
            {"id": "ENST00000000005", "version": 1, "name": "SYNC-201", "biotype": "lncRNA",
             "tags": ["basic", "Ensembl_canonical"],
             "exons": [(12001, 12400), (13001, 13300)], "cds": None,
             "refseq_id": "NR_900005.1"},
        ],
    },
    {
        "gene_id": "ENSG00000000004", "gene_version": 2, "name": "SYND", "strand": "+",
        "biotype": "protein_coding", "ncbi_gene_id": 900004, "description": "synthetic gene D",
        "transcripts": [
            {"id": "ENST00000000006", "version": 3, "name": "SYND-201", "biotype": "protein_coding",
             "tags": ["basic", "Ensembl_canonical"],
             "exons": [(14501, 15600)], "cds": (14601, 15401),
             "protein_id": "ENSP00000000006", "protein_version": 3,
             "refseq_id": "NM_900006.3", "refseq_protein": "NP_900006.2"},
        ],
    },
]

# UniProt entries: (file, header, protein ID of the transcript it matches, or a modified sequence).
# SYNB-201 and SYNB-202 match one SwissProt isoform; SYND has no UniProt entry.
UNIPROT_ENTRIES = [
    ("uniprot.fasta", "sp|P0SYA1|SYNA_HUMAN Synthetic protein A OS=Homo sapiens OX=9606 GN=SYNA PE=1 SV=1", "ENSP00000000001"),
    ("uniprot.fasta", "sp|P0SYB1|SYNB_HUMAN Synthetic protein B OS=Homo sapiens OX=9606 GN=SYNB PE=1 SV=2", "canonical_SYNB"),
    ("uniprot.fasta", "tr|A0A0S1YNA6|A0A0S1YNA6_HUMAN Synthetic protein A OS=Homo sapiens OX=9606 GN=SYNA PE=4 SV=1", "ENSP00000000002"),
    ("uniprot_additional.fasta", "sp|P0SYB1-2|SYNB_HUMAN Isoform 2 of Synthetic protein B OS=Homo sapiens OX=9606 GN=SYNB", "ENSP00000000003"),
]


### Sequence helpers ###

def reverse_complement(seq):
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def translate(seq):
    return "".join(CODON_TABLE[seq[i:i + 3]] for i in range(0, len(seq) - 2, 3))


def wrap(seq, width):
    return "".join(seq[i:i + width] + "\n" for i in range(0, len(seq), width))


def transcript_positions(exons, strand):
    """Genomic positions of the exons in transcript order (5' to 3')."""
    positions = [p for start, end in exons for p in range(start, end + 1)]
    return positions if strand == "+" else positions[::-1]


def cds_positions(transcript, strand):
    """Genomic positions of the CDS and of the stop codon, in transcript order."""
    positions = transcript_positions(transcript["exons"], strand)
    low, high = transcript["cds"]
    first = positions.index(low if strand == "+" else high)
    last = positions.index(high if strand == "+" else low)
    if (last - first + 1) % 3:
        raise ValueError(f"CDS of {transcript['id']} is not a multiple of 3")
    return positions[first:last + 1], positions[last + 1:last + 4]


def segments(positions):
    """Contiguous (start, end) blocks of a list of genomic positions, in the order given."""
    blocks = []
    for p in positions:
        if blocks and abs(p - blocks[-1][-1]) == 1:
            blocks[-1].append(p)
        else:
            blocks.append([p])
    return [(min(b), max(b)) for b in blocks]


def build_genome(rng):
    """
    Random background, canonical splice sites for every intron and an open reading frame for
    every CDS. Bases shared by several transcripts are written once and then kept, so the
    transcripts of one gene stay consistent.
    """
    # Index 0 is padding so that string positions are the 1-based genome coordinates.
    genome = ["N"] + [rng.choice("ACGT") for _ in range(GENOME_LENGTH)]
    locked = {}

    def set_base(pos, base, strand):
        base = base if strand == "+" else COMPLEMENT[base]
        if locked.get(pos, base) != base:
            raise ValueError(f"Conflicting base at position {pos}")
        genome[pos] = locked[pos] = base

    for gene in GENES:
        # GT..AG introns read CT..AC in genome order on the - strand.
        left, right = ("GT", "AG") if gene["strand"] == "+" else ("CT", "AC")
        for transcript in gene["transcripts"]:
            exons = transcript["exons"]
            for (_, prev_end), (next_start, _) in zip(exons[:-1], exons[1:]):
                for offset in range(2):
                    set_base(prev_end + 1 + offset, left[offset], "+")
                    set_base(next_start - 2 + offset, right[offset], "+")

    for gene in GENES:
        strand = gene["strand"]
        for transcript in gene["transcripts"]:
            if not transcript["cds"]:
                continue
            cds, stop = cds_positions(transcript, strand)
            orf = cds + stop
            for i in range(0, len(orf), 3):
                codon_pos = orf[i:i + 3]
                options = ["ATG"] if i == 0 else STOP_CODONS if i == len(orf) - 3 else SENSE_CODONS
                fixed = {k: locked[p] if strand == "+" else COMPLEMENT[locked[p]]
                         for k, p in enumerate(codon_pos) if p in locked}
                options = [c for c in options if all(c[k] == b for k, b in fixed.items())]
                if not options:
                    raise ValueError(f"No valid codon at {codon_pos} for {transcript['id']}")
                codon = rng.choice(options)
                for k, p in enumerate(codon_pos):
                    set_base(p, codon[k], strand)
    return "".join(genome)


def spliced(genome, positions, strand):
    seq = "".join(genome[p] for p in positions)
    return seq if strand == "+" else "".join(COMPLEMENT[b] for b in seq)


### GTF writers ###

def gtf_line(chrom, source, feature, start, end, strand, frame, attributes):
    attr = " ".join(f'{key} "{value}";' for key, value in attributes)
    return f"{chrom}\t{source}\t{feature}\t{start}\t{end}\t.\t{strand}\t{frame}\t{attr}\n"


def transcript_features(transcript, strand):
    """
    Exon, CDS, start/stop codon and UTR blocks of a transcript, in transcript order, as
    (feature, start, end, frame, exon_number).
    """
    exons = sorted(transcript["exons"], reverse=strand == "-")
    features = [("exon", s, e, ".", n) for n, (s, e) in enumerate(exons, 1)]
    utr5, utr3 = [], []
    if transcript["cds"]:
        cds, stop = cds_positions(transcript, strand)
        exon_of = {p: n for n, (s, e) in enumerate(exons, 1) for p in range(s, e + 1)}
        done = 0
        for s, e in segments(cds):
            features.append(("CDS", s, e, (3 - done % 3) % 3, exon_of[s]))
            done += e - s + 1
        start = segments(cds[:3])[0]
        features.append(("start_codon", start[0], start[1], 0, exon_of[cds[0]]))
        stop_block = segments(stop)[0]
        features.append(("stop_codon", stop_block[0], stop_block[1], 0, exon_of[stop[0]]))
        coding = set(cds + stop)
        positions = transcript_positions(transcript["exons"], strand)
        first = positions.index(cds[0])
        utr5 = [p for p in positions[:first] if p not in coding]
        utr3 = [p for p in positions[first:] if p not in coding]
    order = {"exon": 0, "CDS": 1, "start_codon": 2, "stop_codon": 3}
    features.sort(key=lambda f: (f[4], order[f[0]]))
    features += [("five_prime_utr", s, e, ".", None) for s, e in segments(utr5)]
    features += [("three_prime_utr", s, e, ".", None) for s, e in segments(utr3)]
    return features


def write_ensembl_gtf(path):
    exon_ids = {}
    with open(path, "w") as out:
        out.write(f"#!genome-build {ASSEMBLY}\n#!genome-version {ASSEMBLY}\n#!genome-date 2026-10\n"
                  f"#!genome-build-accession GCA_000000000.0\n#!genebuild-last-updated 2026-10\n")
        for gene in GENES:
            strand = gene["strand"]
            source = "ensembl_havana" if gene["biotype"] == "protein_coding" else "havana"
            starts = [s for t in gene["transcripts"] for s, _ in t["exons"]]
            ends = [e for t in gene["transcripts"] for _, e in t["exons"]]
            gene_attr = [("gene_id", gene["gene_id"]), ("gene_version", gene["gene_version"]),
                         ("gene_name", gene["name"]), ("gene_source", source), ("gene_biotype", gene["biotype"])]
            out.write(gtf_line(CHROM, source, "gene", min(starts), max(ends), strand, ".", gene_attr))
            for t in gene["transcripts"]:
                head = gene_attr[:2] + [("transcript_id", t["id"]), ("transcript_version", t["version"])]
                tail = gene_attr[2:] + [("transcript_name", t["name"]), ("transcript_source", source),
                                        ("transcript_biotype", t["biotype"])] + [("tag", tag) for tag in t["tags"]]
                out.write(gtf_line(CHROM, source, "transcript", t["exons"][0][0], t["exons"][-1][1], strand, ".",
                                   head + tail + [("transcript_support_level", "1")]))
                for feature, start, end, frame, exon_number in transcript_features(t, strand):
                    attr = head + ([("exon_number", exon_number)] if exon_number else []) + tail
                    if feature == "exon":
                        exon_id = exon_ids.setdefault((start, end, strand), f"ENSE{len(exon_ids) + 1:011d}")
                        attr += [("exon_id", exon_id), ("exon_version", 1)]
                    elif feature == "CDS":
                        attr += [("protein_id", t["protein_id"]), ("protein_version", t["protein_version"])]
                    out.write(gtf_line(CHROM, source, feature, start, end, strand, frame, attr))


def write_refseq_gtf(path):
    """RefSeq GTF as prepare_refseq_gtf leaves it: only NC_ lines, no header."""
    with open(path, "w") as out:
        for gene in GENES:
            strand = gene["strand"]
            coding = gene["biotype"] == "protein_coding"
            starts = [s for t in gene["transcripts"] for s, _ in t["exons"]]
            ends = [e for t in gene["transcripts"] for _, e in t["exons"]]
            xref = [("db_xref", f"GeneID:{gene['ncbi_gene_id']}")]
            out.write(gtf_line(REFSEQ_CHROM, "BestRefSeq", "gene", min(starts), max(ends), strand, ".",
                               [("gene_id", gene["name"]), ("transcript_id", "")] + xref +
                               [("description", gene["description"]), ("gbkey", "Gene"), ("gene", gene["name"]),
                                ("gene_biotype", gene["biotype"])]))
            for n, t in enumerate(gene["transcripts"], 1):
                head = [("gene_id", gene["name"]), ("transcript_id", t["refseq_id"])] + xref
                biotype = "mRNA" if coding else "lnc_RNA"
                product = f"{gene['description']}, transcript variant {n}" if coding else gene["description"]
                out.write(gtf_line(REFSEQ_CHROM, "BestRefSeq", "transcript", t["exons"][0][0], t["exons"][-1][1], strand, ".",
                                   head + [("gbkey", "mRNA" if coding else "ncRNA"), ("gene", gene["name"]),
                                           ("product", product), ("transcript_biotype", biotype)]))
                features = [f for f in transcript_features(t, strand) if not f[0].endswith("_utr")]
                for feature, start, end, frame, exon_number in sorted(features, key=lambda f: f[0] != "exon"):
                    if feature == "exon":
                        attr = head + [("gene", gene["name"]), ("product", product), ("transcript_biotype", biotype)]
                    else:
                        attr = head + [("gbkey", "CDS"), ("gene", gene["name"]),
                                       ("product", f"synthetic protein {gene['name'][-1]} isoform {chr(96 + n)}"),
                                       ("protein_id", t["refseq_protein"])]
                    out.write(gtf_line(REFSEQ_CHROM, "BestRefSeq", feature, start, end, strand, frame,
                                       attr + [("exon_number", exon_number)]))


### FASTA writers ###

def write_fastas(outdir, genome):
    proteins = {}
    with open(os.path.join(outdir, "genome.fa"), "w") as out:
        out.write(f">{CHROM} dna:chromosome chromosome:{ASSEMBLY}:{CHROM}:1:{GENOME_LENGTH}:1 REF\n")
        out.write(wrap(genome[1:], 60))

    with open(os.path.join(outdir, "ensembl_cdna.fa"), "w") as cdna, \
         open(os.path.join(outdir, "ensembl_pep.fa"), "w") as pep, \
         open(os.path.join(outdir, "refseq_rna.fna"), "w") as rna, \
         open(os.path.join(outdir, "refseq_protein.faa"), "w") as faa:
        for gene in GENES:
            strand = gene["strand"]
            strand_num = 1 if strand == "+" else -1
            gene_desc = f"gene:{gene['gene_id']}.{gene['gene_version']} gene_biotype:{gene['biotype']}"
            symbol_desc = f"gene_symbol:{gene['name']} description:{gene['description']} [Source:Synthetic;Acc:SYN:{gene['ncbi_gene_id']}]"
            for n, t in enumerate(gene["transcripts"], 1):
                exons = t["exons"]
                seq = spliced(genome, transcript_positions(exons, strand), strand)
                location = f"chromosome:{ASSEMBLY}:{CHROM}:{exons[0][0]}:{exons[-1][1]}:{strand_num}"
                cdna.write(f">{t['id']}.{t['version']} cdna {location} {gene_desc} transcript_biotype:{t['biotype']} {symbol_desc}\n")
                cdna.write(wrap(seq, 60))
                if t["cds"]:
                    rna.write(f">{t['refseq_id']} Homo sapiens {gene['description']} ({gene['name']}), transcript variant {n}, mRNA\n")
                else:
                    rna.write(f">{t['refseq_id']} Homo sapiens {gene['description']} ({gene['name']}), long non-coding RNA\n")
                rna.write(wrap(seq, 80))
                if not t["cds"]:
                    continue
                cds, _ = cds_positions(t, strand)
                protein = translate(spliced(genome, cds, strand))
                if protein[0] != "M" or "*" in protein:
                    raise ValueError(f"Bad ORF for {t['id']}")
                proteins[t["protein_id"]] = protein
                cds_location = f"chromosome:{ASSEMBLY}:{CHROM}:{t['cds'][0]}:{t['cds'][1]}:{strand_num}"
                pep.write(f">{t['protein_id']}.{t['protein_version']} pep {cds_location} {gene_desc} "
                          f"transcript:{t['id']}.{t['version']} transcript_biotype:{t['biotype']} {symbol_desc}\n")
                pep.write(wrap(protein, 60))
                faa.write(f">{t['refseq_protein']} synthetic protein {gene['name'][-1]} isoform {chr(96 + n)} [Homo sapiens]\n")
                faa.write(wrap(protein, 80))

    # Canonical SwissProt entry of SYNB: the isoform sequence with 10 extra residues, so it matches no transcript.
    synb = proteins["ENSP00000000003"]
    proteins["canonical_SYNB"] = synb[:40] + "GSGSGSGSGS" + synb[40:]
    for name in ("uniprot.fasta", "uniprot_additional.fasta"):
        with open(os.path.join(outdir, name), "w") as out:
            for file_name, header, source in UNIPROT_ENTRIES:
                if file_name == name:
                    out.write(f">{header}\n{wrap(proteins[source], 60)}")


def main():
    parser = argparse.ArgumentParser(description="Write the synthetic IsoAnnot test dataset.")
    parser.add_argument("--outdir", default=FIXTURE_DIR, help="Output folder (default: %(default)s)")
    args = parser.parse_args()
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    genome = build_genome(random.Random(SEED))
    write_fastas(outdir, genome)
    write_ensembl_gtf(os.path.join(outdir, "ensembl.gtf"))
    write_refseq_gtf(os.path.join(outdir, "refseq.gtf"))
    with open(os.path.join(outdir, "chr2acc"), "w") as out:
        out.write(f"#Chromosome\tAccession.version\n{CHROM}\t{REFSEQ_CHROM}\n")

    # Derived files, made by the pipeline scripts as in rules run_refsqanti and parse_protein_databases.
    def path(name):
        return os.path.join(outdir, name)

    subprocess.run([sys.executable, os.path.join(SCRIPTS_DIR, "referenceSQANTI.py"),
                    "--gtf_file", path("ensembl.gtf"), "--reference_file", path("genome.fa"),
                    "--chr_ref", path("chr2acc"), "--database", "ensembl",
                    "--output_classification", path("ensembl_sqanti_classification.txt"),
                    "--output_junctions", path("ensembl_sqanti_junctions.txt"),
                    "--output_nmd", path("ensembl_sqanti_NMD.txt")], cwd=SCRIPTS_DIR, check=True)
    subprocess.run([sys.executable, os.path.join(SCRIPTS_DIR, "parse_protein_databases.py"),
                    "--uniprot_fasta", path("uniprot.fasta"), path("uniprot_additional.fasta"),
                    "--ensembl_fasta", path("ensembl_pep.fa"), "--refseq_fasta", path("refseq_protein.faa"),
                    "--output", path("parsed_databases.json")], cwd=SCRIPTS_DIR, check=True, stdout=subprocess.DEVNULL)


if __name__ == "__main__":
    main()
