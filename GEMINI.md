# IsoAnnot - Project Context

## Overview
IsoAnnot is a Snakemake-based bioinformatics pipeline designed to generate functional and structural annotations at the isoform level. It integrates data from genomic and proteomic databases (such as Ensembl, RefSeq, UniProt, and Reactome) to categorize and describe transcripts, mapping features for both transcript and protein levels.

## Pipeline Architecture & Core Tools
- **Orchestration**: Snakemake
- **Execution Script**: [isoannot.sh](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/isoannot.sh)
- **Environments**: Conda (`envs/`)
- **Key Modules**:
  - **SQANTI**: Transcripts categorization and cDNA correction mapping.
  - **UtrScan**: Identification of UTR motifs.
  - **RepeatMasker**: Annotation of repetitive elements.
  - **InterProScan**: Functional domain predictions.
  - **UniProt / PhosphoSitePlus**: Incorporation of post-translational modification sites and functional domains.

## Directory Map
- [IsoAnnot/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/): Primary source code directory.
  - [config/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/config/): Contains Snakemake workflow configuration files.
    - [ensembl/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/config/ensembl/): Ensembl reference configurations.
    - [refseq/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/config/refseq/): RefSeq reference configurations.
    - [mytranscripts/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/config/mytranscripts/): Configurations for custom transcript annotations.
    - [generic/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/config/generic/): Unified Snakemake scripts (`Snakefile.smk`, `Snakefile_mytranscripts.smk`, etc.).
  - [envs/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/envs/): Conda environment definition YAML files.
  - [scripts/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/scripts/): Internal Python and shell processing scripts.
  - [software/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/software/): Locally bundled or linked external software (e.g. UtrScan, InterProScan).
  - [data/](file:///home/pabloati/Programs/IsoAnnot/IsoAnnot/data/): Run outputs, temporary files, and downloaded databases.
- [README.md](file:///home/pabloati/Programs/IsoAnnot/README.md): Official documentation.

## Setup & Execution
1. Install Conda / Mamba.
2. Initialize the Snakemake environment:
   ```bash
   conda env create -f IsoAnnot/snakemake.yaml
   conda activate snakemake
   ```
3. Basic Run Command:
   ```bash
   ./isoannot.sh --database <database_option> --species <species_name> [--config option1=value ...]
   ```
   *Note: IsoAnnot is explicitly designed for single-node execution (e.g. within a single SLURM job allocation) using local multi-threading managed by Snakemake (`isoannot.sh`).*
