import glob
import os
import sys

import pytest
import yaml

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, main_path)

CONFIG_DIR = os.path.join(os.path.dirname(main_path), "config")
REFSEQ_KEYS = ["refseq_proteins", "refseq_cdna", "refseq_gtf", "refseq_chr_accessions"]
SPECIES_IN_BOTH = sorted(set(os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(CONFIG_DIR, "ensembl", "*", "config.yaml")))
                         & set(os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(CONFIG_DIR, "refseq", "*", "config.yaml"))))


def read_config(database, species):
    with open(os.path.join(CONFIG_DIR, database, species, "config.yaml")) as h:
        return yaml.safe_load(h)


def assembly_folder(url):
    """The NCBI assembly folder of a RefSeq URL, e.g. GCF_000002035.6_GRCz11."""
    return next(part for part in url.split("/") if part.startswith("GCF_") and "." in part and not part.endswith(".gz"))


@pytest.mark.parametrize("species", SPECIES_IN_BOTH)
def test_ensembl_and_refseq_configs_use_same_assemblies(species):
    """The genome is always the Ensembl one, so refseq mode must use the RefSeq files of the same assembly as ensembl mode."""
    ensembl, refseq = read_config("ensembl", species), read_config("refseq", species)
    for key in ["ensembl_reference"] + REFSEQ_KEYS:
        assert refseq[key] == ensembl[key], key


def test_drerio_refseq_is_grcz11():
    """config/refseq/drerio used the GRCz12tu RefSeq annotation with the GRCz11 Ensembl genome."""
    config = read_config("refseq", "drerio")
    assert "GRCz11" in os.path.basename(config["ensembl_reference"])
    assert {assembly_folder(config[key]) for key in REFSEQ_KEYS} == {"GCF_000002035.6_GRCz11"}
