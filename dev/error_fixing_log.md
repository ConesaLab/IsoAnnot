# IsoAnnot - Protocol Modernization & Critical Error Fix Log

**Date**: July 28, 2026  
**Target Directory**: [IsoAnnot/config/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/config/) and [IsoAnnot/scripts/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/scripts/)

---

## 1. Executive Summary

This document logs the rationale, diagnosis, and refactoring executed to resolve critical protocol and network configuration errors in IsoAnnot. 

The primary fix resolves pipeline crashes on High-Performance Computing (HPC) clusters and cloud environments caused by legacy `ftp://` download URLs, malformed protocol prefixes, and unencrypted `http://` API endpoints.

---

## 2. Diagnosis & Problem Statement

### A. The Legacy `ftp://` Problem on HPC Clusters
* **Root Cause**: Major genomic databases (Ensembl, NCBI RefSeq, UniProt) named their bulk download hostnames with `ftp.` (e.g., `ftp.ensembl.org`, `ftp.ncbi.nlm.nih.gov`, `ftp.uniprot.org`) in the 1990s/2000s. Early configuration files retained `ftp://` schemes.
* **Failure Mechanism**:
  1. **Port 21 & High Port Blocks**: FTP relies on **TCP Port 21** for command signals and a dynamic range of high ports for passive data transfers. HPC compute nodes (SLURM, SGE) and cloud infrastructure block outbound Port 21 and dynamic ports by default for security.
  2. **Proxy Incompatibility**: Enterprise and institutional HTTP proxies (`http_proxy` / `https_proxy`) do not route raw FTP control traffic.
  3. **Snakemake Rule Hangs**: When Snakemake invokes `wget` or `curl` on an `ftp://` URL inside a firewalled node, the command hangs indefinitely until socket timeout or crashes with `Connection refused` / `Network unreachable`.

### B. Malformed Protocol Prefixes
* **Root Cause**: Manual editing attempts in `IsoAnnot/config/ensembl/drerio/config.yaml` resulted in corrupted URL strings such as:
  `ftp://https://ftp.ncbi.nlm.nih.gov/genomes/...`
* **Failure Mechanism**: `wget` and `curl` fail immediately to resolve the hostname `https:`.

### C. Insecure BioMart & External HTTP Endpoints
* **Root Cause**: Python scripts (`layer_go.py`, `layer_reactome.py`, `uniprotPhosphosite_genomicCoordinates.py`, `uniprotPhosphosite_annotation.py`) had hardcoded fallback CLI arguments pointing to unencrypted `http://www.ensembl.org`.
* **Failure Mechanism**: Modern Ensembl web servers enforce 301 Permanent Redirects to `https://www.ensembl.org`. Direct HTTP calls without SSL negotiation can break certain Python REST clients or trigger HTTP 301 unhandled exceptions in strict network environments.

---

## 3. Standardized Transformations Applied

| Original Pattern / URL | Standardized Target | Rationale |
| :--- | :--- | :--- |
| `ftp://ftp.ensembl.org/` | `https://ftp.ensembl.org/` | Modern HTTPS web mirror (Port 443) |
| `ftp://ftp.ensemblgenomes.org/` | `https://ftp.ensemblgenomes.org/` | Modern HTTPS web mirror (Port 443) |
| `ftp://ftp.ncbi.nlm.nih.gov/` | `https://ftp.ncbi.nlm.nih.gov/` | Modern HTTPS web mirror (Port 443) |
| `ftp://ftp.uniprot.org/` | `https://ftp.uniprot.org/` | Modern HTTPS web mirror (Port 443) |
| `ftp://https://ftp.ncbi.nlm.nih.gov/` | `https://ftp.ncbi.nlm.nih.gov/` | Repair corrupted protocol string |
| `http://mirwalk.umm.uni-heidelberg.de/` | `https://mirwalk.umm.uni-heidelberg.de/` | Upgrade to encrypted HTTPS |
| `http://www.ensembl.org` | `https://www.ensembl.org` | Upgrade BioMart CLI default host |

---

## 4. Summary of Refactored Files

### Species Configurations (`IsoAnnot/config/`)
- `config/ensembl/athaliana/config.yaml`
- `config/ensembl/dmelanogaster/config.yaml`
- `config/ensembl/drerio/config.yaml`
- `config/ensembl/hsapiens/config.yaml`
- `config/ensembl/mmusculus/config.yaml`
- `config/refseq/athaliana/config.yaml`
- `config/refseq/dmelanogaster/config.yaml`
- `config/refseq/drerio/config.yaml`
- `config/refseq/hsapiens/config.yaml`
- `config/refseq/mmusculus/config.yaml`
- `config/mytranscripts/hsapiens/config.yaml`
- `config/generic/config.yaml`

### Python Processing & Layer Scripts (`IsoAnnot/scripts/`)
- `scripts/IsoAnnot.py` (Core `ChromosomeMap` module)
- `scripts/referenceSQANTI.py`
- `scripts/mirna_bs_genomic_coord.py`
- `scripts/get_mirna_bs_annotation.py`
- `scripts/uniprotPhosphosite_genomicCoordinates.py`
- `scripts/uniprotPhosphosite_annotation.py`
- `scripts/transcript2reference.py`
- `scripts/layer_exons.py`
- `scripts/layer_go.py`
- `scripts/layer_reactome.py`

---

## 5. Verification Protocol

To verify that these changes function cleanly:
1. **YAML Validation**: Parse all `.yaml` files with `pyyaml` to confirm no syntax or indentation errors were introduced.
2. **Python Compilation**: Execute `python3 -m py_compile` on modified scripts.
3. **HTTP Head Checks**: Execute `curl -I` on sample replaced HTTPS URLs to confirm HTTP 200 OK server response.

---

## 6. Issue 1.3: Chromosome Accession Mapping & Silent Data Loss Prevention

### A. Rationale & Architecture
* **Problem**: In multi-database annotation pipelines, reference files and external feature databases use differing chromosome nomenclature (`1` vs `chr1` vs `NC_000001.11`). Naive single-direction dictionary lookups or version stripping caused unmapped contigs to fail silently, leading to empty output tables.
* **Solution**: Introduced `ChromosomeMap` in [`IsoAnnot/scripts/IsoAnnot.py`](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/scripts/IsoAnnot.py).

### B. Tiered Lookup Hierarchy
1. **Tier 1 (Exact Match)**: Direct key lookup. Prevents collisions between custom/draft scaffolds (`scaffold_1.1` vs `scaffold_1.2`).
2. **Tier 2 (Case & Prefix Normalization)**: Handles capitalization and optional `chr` prefixes (`Chr1`, `chr1`, `CHR1` $\rightarrow$ `1`).
3. **Tier 3 (Accession Version Stripping)**: Version dots (`NC_000001.11` $\rightarrow$ `NC_000001`) are stripped **only as a final fallback** for official NCBI accessions (`NC_`, `NW_`, `NT_`, `AC_`, `NZ_`).

### C. Leading Database Modes
* **Default (`leading_db="ensembl"`)**: Canonical output chromosome target is Ensembl identifiers (`1`).
* **RefSeq Mode (`leading_db="refseq"`)**: Canonical output chromosome target is RefSeq accessions (`NC_000001.11`).

### D. Safeguards & Unit Verification
* **Unmapped Tracking**: `ChromosomeMap` tracks all unmapped queries and total query counts.
* **Safety Threshold**: Scripts check `len(chr_map.unmapped_log) / chr_map.total_queries > 0.50` and raise `RuntimeError` to halt execution if over 50% of feature lookups fail.
* **Unit Tests**: Added [`IsoAnnot/tests/test_chromosome_map.py`](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/tests/test_chromosome_map.py) covering all tiers, leading DB modes, custom scaffolds, and fallback logging.

---

## 7. Issue 1.2: BioMart Resilience (5-Attempt Exponential Backoff Retry & Layer Deactivation)

### A. Rationale & Architecture
* **Problem**: Transient network issues (`HTTP 502 / 504 Gateway Timeout`) or Ensembl server load spikes caused single-attempt `pybiomart` queries to abort long-running pipeline runs.
* **Solution**: Implemented `query_biomart_with_retry` in [`IsoAnnot/scripts/IsoAnnot.py`](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/scripts/IsoAnnot.py#L149-L170).

### B. Retry Mechanism
* **Attempt Count**: 5 maximum attempts.
* **Backoff Strategy**: Exponential delays (5s, 10s, 20s, 40s) between retries.
* **Fatal Failure Behavior**: If the 5th attempt fails, raises a descriptive `RuntimeError` containing explicit instructions for the user on how to deactivate the layer if BioMart remains unreachable.

### C. Layer Deactivation Options for Users
If BioMart is down or network access is restricted, users can bypass BioMart-dependent layers (`layer_go`, `layer_reactome`) without editing code:
1. **Command Line (`isoannot.sh`)**:
   ```bash
   ./isoannot.sh --database ensembl --species hsapiens --config layer_go=no layer_reactome=no
   ```
2. **Species Config (`config.yaml`)**:
   ```yaml
   layer_go: "no"
   layer_reactome: "no"
   ```

### D. Verification
* Unit tests in [`IsoAnnot/tests/test_biomart_retry.py`](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/tests/test_biomart_retry.py) verify successful queries on retry, max attempt limits, and fatal exception messages (**3/3 tests passing**).


