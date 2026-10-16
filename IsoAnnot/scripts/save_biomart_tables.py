#!/usr/bin/env python3
'''
Saves the BioMart answers that layer_go and layer_reactome need (gene ID, gene name and GO terms)
in one table. Set the biomart_table config key to that file and both layers read it instead of
querying BioMart: runs then work without internet and do not change when Ensembl updates BioMart.
Use the archive host of the Ensembl release in your config, e.g. sep2025.archive.ensembl.org for
release 115, since www.ensembl.org always serves the newest release.
'''

import argparse, sys, os, io, logging
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


def https_dataset(host, dataset, schema="default"):
    """
    Dataset-like object whose query() asks BioMart over https. pybiomart 0.2.0 always connects with
    http://, and the Ensembl archive hosts refuse plain http (403), so a release archive can't be
    reached through it. Without a scheme, the host is taken as https.
    """
    url = (host if host.startswith(("http://", "https://")) else f"https://{host}").rstrip("/") + "/biomart/martservice"

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
    parser.add_argument('--biomart_host', required=True, help='BioMart host, e.g. sep2025.archive.ensembl.org (https is used)')
    parser.add_argument('--biomart_dataset', required=True, help='BioMart dataset, e.g. hsapiens_gene_ensembl')
    parser.add_argument('--output', required=True, help='Table to write')
    args = parser.parse_args()

    biomart_schema = "plants_mart" if "plants" in args.biomart_host else "default"
    try:
        rows = save_biomart_table(https_dataset(args.biomart_host, args.biomart_dataset, biomart_schema), args.output)
    except Exception as ex:
        logging.error(str(ex), exc_info=True)
        sys.exit(os.EX_SOFTWARE)
    logging.info(f"Wrote {rows} rows from {args.biomart_host} ({args.biomart_dataset}) to {args.output}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
