#!/usr/bin/env python3
'''
Saves the BioMart answers that layer_go and layer_reactome need (gene ID, gene name and GO terms)
in one table. The get_biomart_table rule runs it once per output folder, and both layers read the
table. It can also be run by hand, to give the table to IsoAnnot with the biomart_table config key:
runs then work without internet.
Use the host of the Ensembl release in your config, e.g. e115.ensembl.org for release 115 or
eg62-plants.ensembl.org for Ensembl Plants 62, since www.ensembl.org always serves the newest release.
'''

import argparse, sys, os, io, re, logging
import urllib.parse
import urllib.request
from types import SimpleNamespace
import pandas as pd
from IsoAnnot import query_biomart_with_retry, BIOMART_TABLE_ATTRIBUTES


def biomart_query_xml(dataset, attributes, schema):
    """BioMart XML query for TSV with header and unique rows, ending with a [success] line when complete."""
    fields = "".join(f'<Attribute name="{a}"/>' for a in attributes)
    return ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
            f'<Query virtualSchemaName="{schema}" formatter="TSV" header="1" uniqueRows="1" completionStamp="1" datasetConfigVersion="0.6">'
            f'<Dataset name="{dataset}" interface="default">{fields}</Dataset></Query>')


def biomart_url(host):
    """
    BioMart address of a host, always over https: the Ensembl archive hosts refuse plain http (403).
    The host can be written with or without http:// or https://.
    """
    return "https://" + re.sub(r"^https?://", "", host.strip()).rstrip("/") + "/biomart/martservice"


def biomart_release(host):
    """
    Release that the BioMart of a host serves, read from the name of its gene database
    (ensembl_mart_115, plants_mart_62).

    :return: the release as an integer, or None when the host does not answer or has no such database.
    """
    try:
        with urllib.request.urlopen(f"{biomart_url(host)}?type=registry", timeout=60) as answer:
            registry = answer.read().decode("utf-8")
    except Exception as err:
        logging.info(f"Could not read the release of {host} ({err})")
        return None
    release = re.search(r'database="\w+_mart_(\d+)"', registry)
    return int(release.group(1)) if release else None


def release_of_url(url):
    """Release in an Ensembl download URL (.../pub/release-115/gtf/...), or None when it has none."""
    release = re.search(r"release-(\d+)", url or "")
    return int(release.group(1)) if release else None


def release_mismatch(biomart, files):
    """
    The classic way to lose annotation without an error: gene names and GO terms of one Ensembl
    release matched against the genes of another.

    :return: the warning to show, or None when the releases agree or one of them is not known.
    """
    if biomart is None or files is None or biomart == files:
        return None
    return (f"BioMart serves Ensembl release {biomart}, but the files of the config are release {files}. Genes that changed "
            f"between the two releases can lose their GO terms and pathways. Set biomart_host to the host of release {files} "
            f"(e{files}.ensembl.org, or eg{files}-plants.ensembl.org for plants), if Ensembl still keeps its archive.")


def https_dataset(host, dataset, schema="default"):
    """
    Dataset-like object whose query() asks BioMart over https, with the host written in any form
    (see biomart_url).
    """
    url = biomart_url(host)

    def query(attributes):
        request = f"{url}?{urllib.parse.urlencode({'query': biomart_query_xml(dataset, attributes, schema)})}"
        with urllib.request.urlopen(request, timeout=900) as answer:
            text = answer.read().decode("utf-8").rstrip("\n")
        # BioMart can stop half-way without an error: only an answer ending in [success] is complete
        if not text.endswith("[success]"):
            raise ValueError(f"BioMart answer from {url} is incomplete (no [success] line at the end)")
        return pd.read_csv(io.StringIO(text[:-len("[success]")]), sep="\t", dtype=str)
    return SimpleNamespace(query=query)


def save_biomart_table(dataset, output):
    """
    Queries the gene ID, gene name and GO attributes and writes them as a tab-separated table with
    BioMart's own column names, through <output>.tmp so an interrupted run leaves no partial table.

    :return: number of rows written.
    """
    table = query_biomart_with_retry(dataset=dataset, attributes=BIOMART_TABLE_ATTRIBUTES,
                                     layer_name="save_biomart_tables", disable_hint=False)
    table.to_csv(f"{output}.tmp", sep="\t", index=False)
    os.replace(f"{output}.tmp", output)
    return len(table)


def main():
    parser = argparse.ArgumentParser(description='Save the BioMart table used by layer_go and layer_reactome')
    parser.add_argument('--biomart_host', required=True, help='BioMart host, e.g. e115.ensembl.org (https is used)')
    parser.add_argument('--biomart_dataset', required=True, help='BioMart dataset, e.g. hsapiens_gene_ensembl')
    parser.add_argument('--output', required=True, help='Table to write')
    parser.add_argument('--ensembl_gtf_url', default=None, help='URL of the Ensembl GTF of the config: its release is compared with the release BioMart serves')
    args = parser.parse_args()

    biomart_schema = "plants_mart" if "plants" in args.biomart_host else "default"
    try:
        rows = save_biomart_table(https_dataset(args.biomart_host, args.biomart_dataset, biomart_schema), args.output)
    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        logging.error("The GO and Reactome layers need this table. To run IsoAnnot without BioMart, give a table saved before "
                      "with 'biomart_table=/abs/path/table.tsv', or turn the two layers off with 'layer_go=no reactome=no'.")
        sys.exit(os.EX_SOFTWARE)
    logging.info(f"Wrote {rows} rows from {args.biomart_host} ({args.biomart_dataset}) to {args.output}")
    if args.ensembl_gtf_url:
        warning = release_mismatch(biomart_release(args.biomart_host), release_of_url(args.ensembl_gtf_url))
        if warning:
            logging.warning(warning)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
