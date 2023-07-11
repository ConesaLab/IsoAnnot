# IsoAnnot

IsoAnnot is a new tool for generating functional and structural annotation at isoform level, capable of collecting and integrating information from different databases to categorize and describe each isoform, including functional and structural information for both transcript and protein.

## Requirements

### Computational requirements
The computational requirements to run IsoAnnot may vary depending on the organism of interest and the size of the transcriptome you want to annotate. 

As reference, to annotate a Human transcriptome of 252205 isoforms, IsoAnnot required 8 cores, 12 GB of RAM, 14 GB of disk space and the execution time was 20 hours. 

The number of cores to be used is specified in the last line of the "isoannot.sh" and can be modified by the user.


### Conda
IsoAnnot has been tested under GNU/Linux and is built in Python 3 using Snakemake in combination with conda environments, allowing for automatic dependency management once these two have been installed.

To install conda you can follow the official instructions: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

Once it is done, it must be followed by snakemake installation via conda: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html or you can activate the snakemake conda environment (snakemake.yaml):

```
conda env create -f snakemake.yaml
conda activate snakemake
```

## Installation

The package is distributed as a compressed file containing the proper structure. In order to use it, you only need to untar it to your desired installation folder, provided the requirements are met.

### External sofware used by IsoAnnot

The workflow uses different programs to generate the information and although many of them will be automatically installed by conda, there are some not present and need to be manually installed:

- **UTRScan**: currently bundled with isoannot.
- **Interproscan**: follow instructions in https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html. By default, IsoAnnot assumes that interproscan.sh is located in the subfolder "software/interproscan/". If you have it installed in another directory, you can always create a symbolic link or modify the parameter "interproscan_path" inside the configuration file "config/generic/config.yaml". 

Please take into account that you will need to activate SignalP and TMHMM https://interproscan-docs.readthedocs.io/en/latest/ActivatingLicensedAnalyses.html. For SignalP, you will also have to modify the signalP configuration file "bin/signalp/4.1/signalp" and specify its absolute path in your system.

- **PhosphoSitePlus** files: IsoAnnot retrieves information regarding post traductional modifications from the PhosphoSitePlus database. 
The necessary files are included in the IsoAnnot compressed file, but you will have to download them manually if you want to update them: https://www.phosphosite.org/staticDownloads


## Usage

IsoAnnot works by gathering information from multiple sources, trying to merge them and generate a final GFF3, provided the required config files are available for the requested species.

The origin of some of those files is what IsoAnnot calls "database", and 3 of them are supported:

- ensembl
- refseq
- mytranscripts

Both ensembl and refseq will generate an annotation reference file, and since the files are publicly available, the workflow will automatically download all of them. The database option "mytranscripts" is designed to use our own transcripts files obtained for instance using Pacbio technology, so this probably will be the option you are looking for.

The basic usage of Isoannot is the following:

```
./isoannot.sh --database <database_option\> --species <species_name\> [--config option1=value option2=value ...]
```

The parameters database and species are mandatory, whereas config is optional with reference databases (ensembl and refseq) but requires at least one working file with "mytranscripts" option. The species code should be written in lowercase including only the first letter of the genus, ie. "Homo sapiens" should be "hsapiens".

As a first example, the generation of Homo sapiens reference annotation for Ensembl would be:

```
./isoannot.sh --database ensembl --species hsapiens
```

Once finished, the final GFF3 output should be located in "data/Hsapiens/". The name of the final output file indicates the species and database parameters used: {species}_tappas_{database}_annotation_file.gff3.

For an example of using our own fasta, we will use "mytranscripts" as database option, and populate the "fasta_cdna" option using the config parameter.

```
./isoannot.sh --database mytranscripts --species stuberosum --config fasta_cdna=/path/to/my/fasta/potato.fasta 
```

## Providing custom configuration files

Snakemake configuration files in IsoAnnot use the YAML file format and are organized on a per-species basis.

To enable a new species analysis for a reference, we need to make sure that all files are available. The easiest way is to use an already defined species as a template, taking Homo sapiens for Ensembl, the config file is located in config/ensembl/hsapiens/config.yaml with the following contents:

```yaml
species:  Homo sapiens
species_name: human

biomart_host: http://www.ensembl.org
biomart_dataset: hsapiens_gene_ensembl

ensembl_cdna: ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
ensembl_proteins: ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
ensembl_gtf:  ftp://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.chr.gtf.gz
ensembl_reference_dir: ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 

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
For a conventional species having all the required files in both Refseq and Ensembl databases, the procedure should be as easy as to change human references and links to the other organism. For instance, to adapt the previous template to rat, we would need to modify the prefix, scientific name, common name and other self-explanatory variables, we would also need to browse the Refseq/Ensembl FTP and copy the address for each required file. For Solanum tuberosum, a preliminary file could be like this:

```yaml
species: Solanum tuberosum
species_name: potato

biomart_host: plants.ensembl.org
biomart_dataset: stuberosum_eg_gene

ensembl_cdna: ftp://ftp.ensemblgenomes.org/pub/release-45/plants/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
ensembl_proteins: ftp://ftp.ensemblgenomes.org/pub/release-45/plants/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
ensembl_gtf: ftp://ftp.ensemblgenomes.org/pub/release-45/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.45.gtf.gz
ensembl_reference: ftp://ftp.ensemblgenomes.org/pub/release-45/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz


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

### Indications

- Set `transcript_versioned` to True if you are executing IsoAnnot in refseq or mytranscripts mode, or if you are annotating a plant (EnsemblPlants).
-  The `ensembl_reference` parameter contains the link to download the ensembl reference genome. If the "*dna.primary_assembly.fa.gz" does not exist, download the "*dna.toplevel.fa.gz"

## Regarding Snakemake

- If you delete or edit any files within the "IsoAnnot/data/" folder, the next time you try to run our pipeline, snkamemake will ask you to `unlock` the folder. To do this you have to edit the last line of the "isoannot.sh" script and add "--unlock", run the pipeline, remove the "--unlock" parameter and run the pipeline again.

- If your execution stopped but some of the steps were completed, next time you run IsoAnnot, it will resume the execution without running again the completed steps. You can do a dry-run to see which steps are going to be executed (add "-n" to the last line of the "isoannot.sh") and you can force IsoAnnot to run all the steps again (add "-F to the last line of the "isoannot.sh").
