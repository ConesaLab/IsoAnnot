import csv
import json
import os
import subprocess
import sys

import pytest
from Bio.Seq import Seq

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

from IsoAnnot import read_chr_ref_acc

ENSEMBL_TO_REFSEQ = {
    "ENST00000000001": "NM_900001.1",
    "ENST00000000002": "NM_900002.1",
    "ENST00000000003": "NM_900003.2",
    "ENST00000000004": "NM_900004.1",
    "ENST00000000005": "NR_900005.1",
    "ENST00000000006": "NM_900006.3",
}

GENERATED_FILES = [
    "genome.fa", "ensembl.gtf", "ensembl_cdna.fa", "ensembl_pep.fa",
    "refseq.gtf", "refseq_rna.fna", "refseq_protein.faa", "chr2acc",
    "uniprot.fasta", "uniprot_additional.fasta",
    "ensembl_sqanti_classification.txt", "ensembl_sqanti_junctions.txt", "ensembl_sqanti_NMD.txt",
    "refseq_sqanti_classification.txt", "refseq_sqanti_junctions.txt", "refseq_sqanti_NMD.txt",
    "parsed_databases.json",
]


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


def get_attribute(attributes, key):
    for field in attributes.split(";"):
        field = field.strip()
        if field.startswith(key + " "):
            return field.split('"')[1]
    return None


def read_gtf_models(path):
    """Exons, CDS and stop codon blocks of each transcript, with versioned transcript and protein IDs."""
    models = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            transcript_id = get_attribute(fields[8], "transcript_id")
            if not transcript_id:
                continue
            version = get_attribute(fields[8], "transcript_version")
            model = models.setdefault(transcript_id, {
                "chrom": fields[0], "strand": fields[6], "exon": [], "CDS": [], "stop_codon": [], "protein": None,
                "versioned_id": f"{transcript_id}.{version}" if version else transcript_id})
            if fields[2] in ("exon", "CDS", "stop_codon"):
                model[fields[2]].append((int(fields[3]), int(fields[4])))
            if fields[2] == "CDS":
                protein_id = get_attribute(fields[8], "protein_id")
                protein_version = get_attribute(fields[8], "protein_version")
                model["protein"] = f"{protein_id}.{protein_version}" if protein_version else protein_id
    return models


def read_chr2acc(path):
    with open(path) as handle:
        return dict(line.rstrip("\n").split("\t") for line in handle if not line.startswith("#"))


def spliced(chrom_seq, blocks, strand):
    seq = "".join(chrom_seq[start - 1:end] for start, end in sorted(blocks))
    return seq if strand == "+" else str(Seq(seq).reverse_complement())


def genome_and_models(synthetic_dir, gtf_name):
    genome = read_fasta(synthetic_dir / "genome.fa")
    to_ensembl = {acc: chrom for chrom, acc in read_chr2acc(synthetic_dir / "chr2acc").items()}
    models = read_gtf_models(synthetic_dir / gtf_name)
    for model in models.values():
        model["seq"] = genome[to_ensembl.get(model["chrom"], model["chrom"])]
    return models


@pytest.fixture(scope="module")
def regenerated_dir(tmp_path_factory):
    outdir = tmp_path_factory.mktemp("synthetic")
    script = os.path.join(os.path.dirname(__file__), "fixtures", "synthetic", "make_fixture.py")
    subprocess.run([sys.executable, script, "--outdir", str(outdir)], check=True)
    return outdir


@pytest.mark.parametrize("name", GENERATED_FILES)
def test_fixture_is_up_to_date(name, synthetic_dir, regenerated_dir):
    """The committed files are what make_fixture.py writes with the current scripts."""
    assert (regenerated_dir / name).read_bytes() == (synthetic_dir / name).read_bytes(), \
        f"{name} is stale: run make_fixture.py and commit the result"


@pytest.mark.parametrize("gtf_name,protein_fasta", [
    ("ensembl.gtf", "ensembl_pep.fa"),
    ("refseq.gtf", "refseq_protein.faa"),
])
def test_proteins_are_translations_of_the_cds(synthetic_dir, gtf_name, protein_fasta):
    """Every protein is the translation of its CDS, which starts with ATG and ends before a stop codon."""
    translations = {}
    for model in genome_and_models(synthetic_dir, gtf_name).values():
        if not model["CDS"]:
            continue
        cds = spliced(model["seq"], model["CDS"], model["strand"])
        stop = spliced(model["seq"], model["stop_codon"], model["strand"])
        assert cds.startswith("ATG")
        assert str(Seq(stop).translate()) == "*"
        translations[model["protein"]] = str(Seq(cds).translate())
    assert translations == read_fasta(synthetic_dir / protein_fasta)


@pytest.mark.parametrize("gtf_name,cdna_fasta", [
    ("ensembl.gtf", "ensembl_cdna.fa"),
    ("refseq.gtf", "refseq_rna.fna"),
])
def test_cdna_is_the_spliced_exons(synthetic_dir, gtf_name, cdna_fasta):
    models = genome_and_models(synthetic_dir, gtf_name)
    cdna = {m["versioned_id"]: spliced(m["seq"], m["exon"], m["strand"]) for m in models.values()}
    assert cdna == read_fasta(synthetic_dir / cdna_fasta)


@pytest.mark.parametrize("gtf_name", ["ensembl.gtf", "refseq.gtf"])
def test_introns_have_canonical_splice_sites(synthetic_dir, gtf_name):
    sites = []
    for model in genome_and_models(synthetic_dir, gtf_name).values():
        exons = sorted(model["exon"])
        for (_, prev_end), (next_start, _) in zip(exons[:-1], exons[1:]):
            intron = spliced(model["seq"], [(prev_end + 1, next_start - 1)], model["strand"])
            sites.append(intron[:2] + intron[-2:])
    assert sites == ["GTAG"] * 13


def test_chr2acc_covers_gtf_chromosomes(synthetic_dir):
    chr2acc = read_chr2acc(synthetic_dir / "chr2acc")
    assert chr2acc == {"1": "NC_000001.11"}
    assert set(read_fasta(synthetic_dir / "genome.fa")) == {"1"}
    assert {m["chrom"] for m in read_gtf_models(synthetic_dir / "ensembl.gtf").values()} == {"1"}
    assert {m["chrom"] for m in read_gtf_models(synthetic_dir / "refseq.gtf").values()} == {"NC_000001.11"}
    chrom_map = read_chr_ref_acc(str(synthetic_dir / "chr2acc"), target_db="ensembl")
    assert chrom_map.get("NC_000001.11") == "1"


def test_refseq_models_match_ensembl(synthetic_dir):
    """The RefSeq GTF has the same transcript models as the Ensembl GTF, on NC_000001.11."""
    ensembl = read_gtf_models(synthetic_dir / "ensembl.gtf")
    refseq = read_gtf_models(synthetic_dir / "refseq.gtf")
    assert sorted(refseq) == sorted(ENSEMBL_TO_REFSEQ.values())
    for ensembl_id, refseq_id in ENSEMBL_TO_REFSEQ.items():
        for key in ("strand", "exon", "CDS", "stop_codon"):
            assert ensembl[ensembl_id][key] == refseq[refseq_id][key]


def test_identical_protein_pair(synthetic_dir):
    """SYNB-201 and SYNB-202 differ in their 5'UTR but encode the same protein."""
    models = read_gtf_models(synthetic_dir / "ensembl.gtf")
    assert models["ENST00000000003"]["exon"] != models["ENST00000000004"]["exon"]
    assert models["ENST00000000003"]["CDS"] == models["ENST00000000004"]["CDS"]
    proteins = read_fasta(synthetic_dir / "ensembl_pep.fa")
    same = [name for name, seq in proteins.items() if seq == proteins["ENSP00000000003.2"]]
    assert same == ["ENSP00000000003.2", "ENSP00000000004.1"]


def test_classification_matches_gtf(synthetic_dir):
    """The ensembl classification has one row per transcript and unversioned refProt, as in a real run."""
    models = read_gtf_models(synthetic_dir / "ensembl.gtf")
    with open(synthetic_dir / "ensembl_sqanti_classification.txt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["isoform"] for row in rows] == sorted(ENSEMBL_TO_REFSEQ)
    for row in rows:
        model = models[row["isoform"]]
        assert (row["chrom"], row["strand"]) == ("1", model["strand"])
        assert int(row["exons"]) == len(model["exon"])
        assert int(row["length"]) == sum(end - start + 1 for start, end in model["exon"])
        if not model["CDS"]:
            assert (row["coding"], row["refProt"], row["CDS_start"]) == ("non_coding", "NA", "NA")
            continue
        assert row["coding"] == "coding"
        assert row["refProt"] == model["protein"].split(".")[0]
        assert int(row["CDS_length"]) == sum(end - start + 1 for start, end in model["CDS"])
    cds_start = {row["isoform"]: int(row["CDS_start"]) for row in rows if row["coding"] == "coding"}
    assert cds_start == {"ENST00000000001": 101, "ENST00000000002": 101, "ENST00000000003": 51,
                         "ENST00000000004": 111, "ENST00000000006": 101}


def test_nmd_and_junction_tables(synthetic_dir):
    with open(synthetic_dir / "ensembl_sqanti_NMD.txt") as handle:
        nmd = dict(line.rstrip("\n").split("\t") for line in handle)
    assert nmd == {"ENST00000000001": "nonNMD", "ENST00000000002": "NMD", "ENST00000000003": "nonNMD",
                   "ENST00000000004": "nonNMD", "ENST00000000005": "NA", "ENST00000000006": "nonNMD"}
    with open(synthetic_dir / "ensembl_sqanti_junctions.txt") as handle:
        junctions = list(csv.DictReader(handle, delimiter="\t"))
    assert len(junctions) == 13
    assert {(j["splice_site"], j["canonical"]) for j in junctions} == {("GTAG", "canonical")}


def test_parsed_databases_json(synthetic_dir):
    """The identical SYNB protein maps to one SwissProt isoform; SYND has no UniProt entry."""
    pep = read_fasta(synthetic_dir / "ensembl_pep.fa")
    uniprot = read_fasta(synthetic_dir / "uniprot.fasta")
    with open(synthetic_dir / "parsed_databases.json") as handle:
        parsed = json.load(handle)
    assert parsed == {
        pep["ENSP00000000001.2"]: {"swissprot": ["P0SYA1"], "refseq": ["NP_900001"], "ensembl": ["ENSP00000000001"]},
        pep["ENSP00000000002.1"]: {"trembl": ["A0A0S1YNA6"], "refseq": ["NP_900002"], "ensembl": ["ENSP00000000002"]},
        pep["ENSP00000000003.2"]: {"swissprot": ["P0SYB1-2"], "refseq": ["NP_900003", "NP_900004"],
                                   "ensembl": ["ENSP00000000003", "ENSP00000000004"]},
        pep["ENSP00000000006.3"]: {"refseq": ["NP_900006"], "ensembl": ["ENSP00000000006"]},
        uniprot["sp|P0SYB1|SYNB_HUMAN"]: {"swissprot": ["P0SYB1"]},
    }


def test_run_script_runs_a_layer_on_the_fixture(synthetic_dir, run_script, tmp_path):
    """The run_script fixture runs a pipeline script on the fixture, with IsoAnnot.py importable."""
    output = tmp_path / "layer_exons.gtf"
    result = run_script("layer_exons.py", "--gtf_file", synthetic_dir / "ensembl.gtf", "--output", output)
    assert result.returncode == 0, result.stderr
    lines = output.read_text().splitlines()
    assert len(lines) == 19
    assert lines[0] == "ENST00000000001\ttappAS\texon\t1001\t1200\t.\t+\t.\tChr=chr1"
