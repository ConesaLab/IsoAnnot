# IsoAnnot

IsoAnnot is a new tool for generating functional and structural annotation at isoform level, capable of collecting and integrating information from different databases to categorize and describe each isoform, including functional and structural information for both transcript and protein.

## Table of Contents

- [Requirements](#requirements)
  - [Computational Requirements](#computational-requirements)
  - [Software Prerequisites](#software-prerequisites)
- [Installation Prerequisites](#installation-prerequisites)
  - [Installing Conda](#installing-conda)
  - [Installing Snakemake](#installing-snakemake)
  - [External Software](#external-software)
- [Installation](#installation)
- [How to Run IsoAnnot](#how-to-run-isoannot)
  - [Basic Usage](#basic-usage)
  - [Example: Ensembl Reference Annotation](#example-ensembl-reference-annotation)
  - [Example: Using Custom Transcripts](#example-using-custom-transcripts)
  - [Command-line Parameters](#command-line-parameters)
- [Configuration Files](#configuration-files)
  - [Where to Find Config Files](#where-to-find-config-files)
  - [How to Modify Config Files](#how-to-modify-config-files)
  - [Generic Configuration](#generic-configuration)
- [Creating a New Config File](#creating-a-new-config-file)
  - [Step-by-Step Guide](#step-by-step-guide)
  - [Config File Template](#config-file-template)
  - [Configuration Parameters Explained](#configuration-parameters-explained)
- [Output](#output)
  - [Output Structure](#output-structure)
  - [Main Output Files](#main-output-files)
  - [Understanding the GFF3 Annotation File](#understanding-the-gff3-annotation-file)
- [Working with Snakemake](#working-with-snakemake)

---

## Requirements

### Computational Requirements

The computational requirements to run IsoAnnot may vary depending on the organism of interest and the size of the transcriptome you want to annotate. 

**Reference benchmark** (Human transcriptome):
- **Transcriptome size**: 252,205 isoforms
- **CPU cores**: 8 cores
- **Memory**: 12 GB RAM
- **Disk space**: 14 GB
- **Execution time**: ~20 hours

The number of cores can be modified by editing the `--cores` parameter in the last line of `IsoAnnot/isoannot.sh` (default is 8 cores).

### Software Prerequisites

IsoAnnot requires the following software to be installed before use:

1. **Operating System**: GNU/Linux (tested and supported)
2. **Python**: Python 3 (managed automatically by conda)
3. **Conda**: For dependency management
4. **Snakemake**: Workflow management system (version 7.x recommended)

---

## Installation Prerequisites

### Installing Conda

Conda is required for automatic dependency management. IsoAnnot has been tested under GNU/Linux and is built using conda environments.

**Installation steps:**

1. Download and install Miniconda or Anaconda following the official instructions:
   - Official Conda installation guide: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

2. Verify conda installation:
   ```bash
   conda --version
   ```

### Installing Snakemake

After installing conda, you need to install Snakemake. IsoAnnot includes a pre-configured environment file for this purpose.

**Option 1: Using the included snakemake.yaml file** (Recommended)

```bash
cd IsoAnnot
conda env create -f snakemake.yaml
conda activate snakemake
```

**Option 2: Manual installation via conda**

Follow the official Snakemake installation guide: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

### External Software

The workflow uses different programs to generate annotations. While many dependencies are automatically installed by conda, some require manual installation:

#### 1. UTRScan
- **Status**: Currently bundled with IsoAnnot
- **Action required**: None (already included)

#### 2. InterProScan
InterProScan is used for protein domain and functional annotation.

**Installation:**
1. Use the provided installation script:
   ```bash
   cd IsoAnnot
   ./InterproScan_install.sh
   ```

2. Activate SignalP and TMHMM by following the installation script instructions.

**Configuration:**
- Default location: `IsoAnnot/software/interproscan/`
- If installed elsewhere, create a symbolic link or modify the `interproscan_path` parameter in `config/generic/config.yaml`

**Important**: Do not move the IsoAnnot directory after installing InterProScan. If you must relocate it, re-run the InterProScan installation script.

#### 3. PhosphoSitePlus Files
IsoAnnot retrieves post-translational modification information from the PhosphoSitePlus database.

- **Status**: Necessary files are included in the IsoAnnot package
- **Updates**: Download manually from https://www.phosphosite.org/staticDownloads if you need updated files

---

## Installation

IsoAnnot is distributed as a compressed file containing the proper directory structure.

**Installation steps:**

1. Extract the package to your desired installation folder:
   ```bash
   tar -xzf IsoAnnot.tar.gz
   cd IsoAnnot
   ```

2. Ensure all prerequisites are installed (see [Installation Prerequisites](#installation-prerequisites))

3. Install external software (see [External Software](#external-software))

4. Activate the snakemake conda environment:
   ```bash
   conda activate snakemake
   ```

You're now ready to run IsoAnnot!

---

## How to Run IsoAnnot

IsoAnnot works by gathering information from multiple sources, merging them, and generating a final GFF3 annotation file.

### Basic Usage

The basic command structure is:

```bash
cd IsoAnnot
./isoannot.sh --database <database_option> --species <species_name> [--config option1=value option2=value ...]
```

**Required parameters:**
- `--database`: Database source (`ensembl`, `refseq`, or `mytranscripts`)
- `--species`: Species code in lowercase (first letter of genus + species, e.g., `hsapiens` for *Homo sapiens*)

**Optional parameter:**
- `--config`: Additional configuration options (required for `mytranscripts`)

**Supported databases:**
- `ensembl` - Ensembl database reference files
- `refseq` - RefSeq (NCBI) database reference files  
- `mytranscripts` - Your own custom transcript files (e.g., from PacBio sequencing)

Both ensembl and refseq will automatically download all required files. The "mytranscripts" option is designed for custom transcript files.

### Example: Ensembl Reference Annotation

To generate *Homo sapiens* reference annotation from Ensembl:

```bash
cd IsoAnnot
./isoannot.sh --database ensembl --species hsapiens
```

**Output location**: `IsoAnnot/data/Hsapiens/human_tappas_ensembl_annotation_file.gff3`

### Example: Using Custom Transcripts

To annotate your own transcripts (e.g., from PacBio sequencing):

```bash
cd IsoAnnot
./isoannot.sh --database mytranscripts --species stuberosum --config fasta_cdna=/path/to/my/fasta/potato.fasta
```

This example annotates potato (*Solanum tuberosum*) transcripts from a custom FASTA file.

### Command-line Parameters

| Parameter | Description | Required | Values |
|-----------|-------------|----------|--------|
| `--database` | Database source | Yes | `ensembl`, `refseq`, `mytranscripts` |
| `--species` | Species identifier | Yes | Lowercase species code (e.g., `hsapiens`, `mmusculus`) |
| `--config` | Override config values | Conditional* | `key=value` pairs |
| `--configfile` | Custom config file path | No | Path to YAML file |
| `--snakefile` | Custom snakefile path | No | Path to Snakefile |

*Required when using `--database mytranscripts` with custom files

**Species code format**: Write species codes in lowercase using only the first letter of the genus followed by the full species name. Examples:
- *Homo sapiens* → `hsapiens`
- *Mus musculus* → `mmusculus`
- *Drosophila melanogaster* → `dmelanogaster`
- *Solanum tuberosum* → `stuberosum`

---

## Configuration Files

Configuration files control how IsoAnnot processes data for each species and database combination. Snakemake configuration files in IsoAnnot use the YAML file format and are organized on a per-species basis.

### Where to Find Config Files

Configuration files are organized in a hierarchical directory structure:

```
IsoAnnot/config/
├── ensembl/
│   ├── hsapiens/
│   │   ├── config.yaml
│   │   └── Snakefile.smk
│   ├── mmusculus/
│   │   ├── config.yaml
│   │   └── Snakefile.smk
│   └── ...
├── refseq/
│   ├── hsapiens/
│   │   ├── config.yaml
│   │   └── Snakefile.smk
│   └── ...
├── mytranscripts/
│   ├── hsapiens/
│   │   ├── config.yaml
│   │   └── Snakefile.smk
│   └── ...
└── generic/
    ├── config.yaml          # Generic settings
    ├── Snakefile.smk        # Main workflow
    ├── Snakefile_ensembl.smk
    ├── Snakefile_refseq.smk
    └── Snakefile_mytranscripts.smk
```

**Path structure**: `config/<database>/<species>/config.yaml`

**Examples**:
- Human Ensembl: `config/ensembl/hsapiens/config.yaml`
- Mouse RefSeq: `config/refseq/mmusculus/config.yaml`
- Custom human transcripts: `config/mytranscripts/hsapiens/config.yaml`

### How to Modify Config Files

To modify an existing configuration:

1. Navigate to the config file:
   ```bash
   cd IsoAnnot/config/<database>/<species>/
   nano config.yaml
   ```

2. Edit parameters as needed (see [Configuration Parameters Explained](#configuration-parameters-explained))

3. Save the file

4. Run IsoAnnot with the updated configuration:
   ```bash
   cd IsoAnnot
   ./isoannot.sh --database <database> --species <species>
   ```

**Common modifications**:
- Update database URLs to newer releases
- Change file paths for custom data
- Adjust species-specific parameters
- Modify the `transcript_versioned` flag

### Generic Configuration

The generic configuration file (`config/generic/config.yaml`) contains global settings used across all species:

```yaml
interproscan_path: "software/interproscan/interproscan.sh"
pfam_clan_url: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
dir_sqanti: "scripts/sqanti3/"
```

**Key parameters**:
- `interproscan_path`: Path to InterProScan executable
- `pfam_clan_url`: URL for Pfam clan database
- `dir_sqanti`: Directory containing SQANTI3 scripts

---

## Creating a New Config File

To enable IsoAnnot for a new species, you need to create configuration files with links to all required data sources.

### Step-by-Step Guide

1. **Choose a template**: Use an existing species config as a starting point
   ```bash
   cd IsoAnnot/config/<database>/
   mkdir <new_species>
   cp hsapiens/config.yaml <new_species>/
   cp hsapiens/Snakefile.smk <new_species>/
   ```

2. **Find required files**: Browse the appropriate database FTP site to locate files for your species:
   - **Ensembl**: https://ftp.ensembl.org/pub/ (or https://ftp.ensemblgenomes.org/ for non-vertebrates)
   - **RefSeq**: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/

3. **Edit the config file**: Update all species-specific parameters (see template below)

4. **Verify the configuration**: Check that all URLs are valid and accessible

5. **Test the configuration**: Run IsoAnnot with the new configuration

### Config File Template

Here's a complete template for creating a new species configuration. The example below shows the config file format for *Homo sapiens* (human) using Ensembl, which is located in `config/ensembl/hsapiens/config.yaml`:

```yaml
species:  Homo sapiens
species_name: human

biomart_host: http://www.ensembl.org
biomart_dataset: hsapiens_gene_ensembl

ensembl_cdna: ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
ensembl_proteins: ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
ensembl_gtf:  ftp://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.chr.gtf.gz
ensembl_reference: ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 

prefix: Hsapiens
db: ensembl
transcript_versioned: False

refseq_protein_dir: ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/
refseq_protein_fasta: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz
refseq_gtf: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
refseq_chr_accessions: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc

uniprot_fasta:
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606_additional.fasta.gz

uniprot_dat:
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.dat.gz
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606_additional.dat.gz

reactome: https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt
```

For a conventional species having all the required files in both RefSeq and Ensembl databases, the procedure should be as easy as changing human references and links to the other organism. You would need to:
1. Modify the prefix, scientific name, common name
2. Browse the RefSeq/Ensembl FTP sites to find the correct URLs for each required file
3. Update the BioMart dataset identifier

Here's another example for *Solanum tuberosum* (potato):

```yaml
species: Solanum tuberosum
species_name: potato

biomart_host: plants.ensembl.org
biomart_dataset: stuberosum_eg_gene

ensembl_cdna: ftp://ftp.ensemblgenomes.org/pub/release-45/plants/fasta/solanum_tuberosum/cdna/Solanum_tuberosum.Assembly.cdna.all.fa.gz
ensembl_proteins: ftp://ftp.ensemblgenomes.org/pub/release-45/plants/fasta/solanum_tuberosum/pep/Solanum_tuberosum.Assembly.pep.all.fa.gz
ensembl_gtf: ftp://ftp.ensemblgenomes.org/pub/release-45/plants/gtf/solanum_tuberosum/Solanum_tuberosum.Assembly.45.gtf.gz
ensembl_reference: ftp://ftp.ensemblgenomes.org/pub/release-45/plants/fasta/solanum_tuberosum/dna/Solanum_tuberosum.Assembly.dna.toplevel.fa.gz

prefix: Stuberosum
db: ensembl
transcript_versioned: True

refseq_protein_dir: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Solanum_tuberosum/latest_assembly_versions/GCF_000226075.1_SolTub_3.0/
refseq_protein_fasta: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Solanum_tuberosum/latest_assembly_versions/GCF_000226075.1_SolTub_3.0/GCF_000226075.1_SolTub_3.0_protein.faa.gz
refseq_chr_accessions: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Solanum_tuberosum/latest_assembly_versions/GCF_000226075.1_SolTub_3.0/GCF_000226075.1_SolTub_3.0_assembly_structure/Primary_Assembly/scaffold_localID2acc
refseq_gtf: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Solanum_tuberosum/latest_assembly_versions/GCF_000226075.1_SolTub_3.0/GCF_000226075.1_SolTub_3.0_genomic.gtf.gz

uniprot_fasta:
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000011115/UP000011115_4113.fasta.gz
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000011115/UP000011115_4113_additional.fasta.gz

uniprot_dat:
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000011115/UP000011115_4113.dat.gz
  - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000011115/UP000011115_4113_additional.dat.gz

reactome: https://plantreactome.gramene.org/download/current/Ensembl2PlantReactome_All_Levels.txt
```

### Configuration Parameters Explained

| Parameter | Description | Example |
|-----------|-------------|---------|
| `species` | Full scientific name | `Homo sapiens` |
| `species_name` | Common name (lowercase) | `human` |
| `biomart_host` | BioMart server URL | `http://www.ensembl.org` |
| `biomart_dataset` | BioMart dataset identifier | `hsapiens_gene_ensembl` |
| `ensembl_cdna` | Ensembl cDNA FASTA file URL | URL to .fa.gz file |
| `ensembl_proteins` | Ensembl protein FASTA file URL | URL to .fa.gz file |
| `ensembl_gtf` | Ensembl GTF annotation file URL | URL to .gtf.gz file |
| `ensembl_reference` | Ensembl reference genome URL | URL to .fa.gz file |
| `prefix` | Output directory prefix (capitalized) | `Hsapiens` |
| `db` | Database type | `ensembl`, `refseq`, or `mytranscripts` |
| `transcript_versioned` | Whether transcripts have version numbers | `True` or `False` |
| `refseq_protein_dir` | RefSeq protein directory URL | FTP directory URL |
| `refseq_protein_fasta` | RefSeq protein FASTA URL | URL to .faa.gz file |
| `refseq_chr_accessions` | RefSeq chromosome accessions file | URL to accession mapping file |
| `refseq_gtf` | RefSeq GTF annotation file URL | URL to .gtf.gz file |
| `uniprot_fasta` | UniProt FASTA file URLs (list) | List of URLs |
| `uniprot_dat` | UniProt DAT file URLs (list) | List of URLs |
| `reactome` | Reactome pathway mapping file URL | URL or pathway-specific URL |
| `layer_go` | Include Gene Ontology layer | `si` (yes) or omit |

**Important notes**:

1. **`transcript_versioned`**: Set to `True` if:
   - Using RefSeq database
   - Using mytranscripts mode
   - Annotating plants (EnsemblPlants)

2. **`ensembl_reference`**: 
   - Try `*dna.primary_assembly.fa.gz` first
   - If not available, use `*dna.toplevel.fa.gz`

3. **Finding UniProt proteome IDs**:
   - Search for your species at https://www.uniprot.org/proteomes
   - Look for "reference proteome" entries
   - Use the proteome ID in the URLs (e.g., `UP000005640` for human)

4. **Reactome URLs**:
   - Animals: `https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt`
   - Plants: `https://plantreactome.gramene.org/download/current/Ensembl2PlantReactome_All_Levels.txt`

---

## Output

### Output Structure

IsoAnnot generates its output in a structured directory hierarchy within the `IsoAnnot/data/` folder:

```
IsoAnnot/data/
└── <Prefix>/                                    # e.g., Hsapiens/
    ├── <species_name>_tappas_<db>_annotation_file.gff3     # Main output
    ├── <species_name>_tappas_<db>_annotation_file.gff3_mod # Modified GFF3
    ├── config/                                  # Downloaded config files
    │   ├── ensembl/
    │   ├── refseq/
    │   └── global/
    ├── output/
    │   └── <db>/                               # Database-specific outputs
    │       ├── layers/                         # Annotation layers
    │       │   ├── go.gtf
    │       │   ├── interpro.gtf
    │       │   ├── reactome.gtf
    │       │   └── ...
    │       ├── transcripts/                    # Transcript files
    │       ├── proteins/                       # Protein sequences
    │       └── ...
    └── tmp/                                    # Temporary processing files
```

**Directory naming**:
- `<Prefix>`: Capitalized species prefix from config (e.g., `Hsapiens`, `Mmusculus`, `Stuberosum`)
- `<species_name>`: Lowercase common name from config (e.g., `human`, `mouse`, `potato`)
- `<db>`: Database used (e.g., `ensembl`, `refseq`, `mytranscripts`)

### Main Output Files

#### Primary Annotation File

**File**: `<species_name>_tappas_<db>_annotation_file.gff3`

This is the main output file containing comprehensive isoform-level annotations.

**Example**: `human_tappas_ensembl_annotation_file.gff3`

**Location**: `IsoAnnot/data/<Prefix>/`

**Content**: GFF3-formatted annotation with:
- Gene and transcript structures
- Protein-coding predictions
- Functional annotations from multiple databases
- Structural features
- Post-translational modifications

#### Modified Annotation File

**File**: `<species_name>_tappas_<db>_annotation_file.gff3_mod`

A modified version of the main GFF3 file optimized for downstream analysis tools.

### Understanding the GFF3 Annotation File

The output GFF3 file integrates information from multiple sources:

**Structural information**:
- Gene and transcript coordinates
- Exon/intron structure
- CDS (coding sequence) regions
- UTR regions (5' and 3')

**Functional annotations** (in attributes column):
- **Gene Ontology (GO)**: Biological process, molecular function, cellular component
- **InterPro**: Protein domains, families, and functional sites
- **Pfam**: Protein family classifications
- **Reactome**: Pathway associations
- **UniProt**: Protein function descriptions

**Post-translational modifications**:
- Phosphorylation sites
- Other PTMs from PhosphoSitePlus

**Example GFF3 attributes**:
```
gene_id=ENSG00000000003;transcript_id=ENST00000000003;GO=GO:0005515,GO:0003824;
InterPro=IPR001478,IPR015421;Reactome=R-HSA-112316;UniProt=P12345
```

**Using the output**:
- Import into genome browsers (IGV, UCSC Genome Browser)
- Use with tappAS for isoform-level functional analysis
- Parse programmatically for custom analyses
- Filter by specific annotation types

---

## Working with Snakemake

IsoAnnot uses Snakemake for workflow management. Here are common Snakemake operations:

### Unlocking the Working Directory

If you delete or edit files in `IsoAnnot/data/`, Snakemake may lock the directory. To unlock:

1. Edit `IsoAnnot/isoannot.sh` and add `--unlock` to the last line:
   ```bash
   exec snakemake ... --unlock
   ```

2. Run the pipeline:
   ```bash
   ./isoannot.sh --database <db> --species <species>
   ```

3. Remove `--unlock` from `isoannot.sh`

4. Run the pipeline normally

### Resuming Interrupted Runs

If execution stops, IsoAnnot will automatically resume from the last completed step on the next run. Completed steps are not re-executed.

### Dry Run (Preview)

To see which steps will be executed without running them:

1. Edit the last line of `IsoAnnot/isoannot.sh` and add `-n`:
   ```bash
   exec snakemake ... -n
   ```

2. Run the pipeline

3. Remove `-n` to execute normally

### Force Re-run

To force IsoAnnot to re-execute all steps:

1. Edit the last line of `IsoAnnot/isoannot.sh` and add `-F`:
   ```bash
   exec snakemake ... -F
   ```

2. Run the pipeline

3. Remove `-F` for subsequent runs

### Adjusting CPU Cores

To change the number of CPU cores used:

1. Edit the last line of `IsoAnnot/isoannot.sh`
2. Modify the `--cores` parameter (default is 8):
   ```bash
   exec snakemake ... --cores 16  # Use 16 cores
   ```

### Common Snakemake Options

Add these to the last line of `isoannot.sh` after the existing parameters:

| Option | Description |
|--------|-------------|
| `-n` or `--dry-run` | Show what would be done without executing |
| `-F` or `--force` | Force re-execution of all steps |
| `--unlock` | Unlock the working directory |
| `--cores N` | Use N CPU cores |
| `-p` | Print shell commands (already included) |
| `--rerun-incomplete` | Re-run incomplete jobs (already included) |

---

## Troubleshooting

**Problem**: "The snakefile or configfile requested do not exist"
- **Solution**: Ensure config files exist for your species at `config/<database>/<species>/`

**Problem**: InterProScan not found
- **Solution**: Run `./InterproScan_install.sh` or verify `interproscan_path` in `config/generic/config.yaml`

**Problem**: Out of memory errors
- **Solution**: Increase available RAM or reduce the number of cores used

**Problem**: Download errors for database files
- **Solution**: Check internet connection and verify URLs in config file are current

**Problem**: Snakemake directory locked
- **Solution**: Use `--unlock` option (see [Unlocking the Working Directory](#unlocking-the-working-directory))

---

## Citation

If you use IsoAnnot in your research, please cite:

[Citation information to be added]

---

## Support

For issues, questions, or contributions:
- **GitHub Issues**: https://github.com/ConesaLab/IsoAnnot/issues
- **Documentation**: This README

---

## License

[License information to be added]
