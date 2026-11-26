#!/usr/bin/env python3
"""
Writes the sequence of each transcript of a GTF, spliced from the genome FASTA. Used in refseq mode instead of
NCBI's rna.fna, whose sequences can differ from the genome (poly(A) tails, insertions): the transcript
positions of UtrScan and RepeatMasker would then not match the exons of the GFF3, which come from the GTF.
"""

import argparse, logging, os, re, sys
from collections import defaultdict
from IsoAnnot import openfile, read_chr_ref_acc

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
LINE_WIDTH = 60
TRANSCRIPT_ID_REGEX = re.compile(r'transcript_id "([^"]+)"')


def read_exons(gtf_file, chrom_map):
    """
    Exons of each transcript, grouped by the genome name of their chromosome.

    :return: {chromosome: {transcript_id: (strand, [(start, end), ...])}}, in GTF order of first appearance.
    """
    transcripts = defaultdict(dict)
    with openfile(gtf_file, "rt") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            match = TRANSCRIPT_ID_REGEX.search(fields[8])
            if not match:
                continue
            chrom = chrom_map.get(fields[0])
            exons = transcripts[chrom].setdefault(match.group(1), (fields[6], []))[1]
            exons.append((int(fields[3]), int(fields[4])))
    return transcripts


def read_genome(genome_file):
    """Yields (name, sequence) for each chromosome, one at a time, so the whole genome is never in memory."""
    name, chunks = None, []
    with openfile(genome_file, "rt") as fasta:
        for line in fasta:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def write_transcripts(transcripts, genome_file, output):
    """
    :return: number of transcripts written and the transcripts left out because their chromosome is not in
    the genome.
    """
    n_written = 0
    with open(f"{output}.tmp", "w") as out:
        for chrom, chrom_seq in read_genome(genome_file):
            for transcript_id, (strand, exons) in transcripts.pop(chrom, {}).items():
                seq = "".join(chrom_seq[start - 1:end] for start, end in sorted(exons))
                if strand == "-":
                    seq = seq.translate(COMPLEMENT)[::-1]
                out.write(f">{transcript_id}\n")
                for i in range(0, len(seq), LINE_WIDTH):
                    out.write(f"{seq[i:i + LINE_WIDTH]}\n")
                n_written += 1
    os.replace(f"{output}.tmp", output)
    missing = [(chrom, transcript_id) for chrom, chrom_transcripts in transcripts.items() for transcript_id in chrom_transcripts]
    return n_written, missing


def main():
    parser = argparse.ArgumentParser(description="Write the transcript sequences of a GTF, spliced from the genome")
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--genome", required=True, help="Genome FASTA (Ensembl chromosome names)")
    parser.add_argument("--chr_ref", required=True, help="NCBI chr2acc table, to translate the GTF chromosome names to the genome ones")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # The GTF has RefSeq names (NC_000001.11) and the genome Ensembl ones (1)
    chrom_map = read_chr_ref_acc(args.chr_ref, target_db="ensembl")
    transcripts = read_exons(args.gtf, chrom_map)
    chrom_map.check_unmapped("transcripts_from_genome.py")
    n_written, missing = write_transcripts(transcripts, args.genome, args.output)
    logging.info(f"{n_written} transcript sequences written to {args.output}")
    if missing:
        chroms = sorted({chrom for chrom, _ in missing})
        logging.warning(f"{len(missing)} transcripts are on chromosomes that are not in the genome ({', '.join(chroms[:10])}) "
                        f"and have no sequence, so UtrScan and RepeatMasker do not annotate them")
    if n_written == 0:
        logging.error(f"No transcript of {args.gtf} is on a chromosome of {args.genome}. Check the chromosome names and {args.chr_ref}")
        sys.exit(1)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
