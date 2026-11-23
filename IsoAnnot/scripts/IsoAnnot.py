"""
IsoAnnot functions
Modified by Alessandra Martinez
"""
import gzip, itertools, logging, re
from collections import Counter, defaultdict
from operator import itemgetter

strand_table = {
    '+': +1,
    '-': -1,
}

def argparse_nullable(arg_value):
    """
    Return arguments when exist.

    Args:
        arg_value: argument value

    Returns:
        arg_value if exists, None if it doesn't
    """
    if not arg_value or not len(arg_value):
        return None
    return arg_value


def openfile(filename, mode='r'):
    """
    Open files compressed (gz) or not.

    Args:
        filename (str): name of file.
        mode (str): opening mode. Read as default.

    Returns:
        open file
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
        #return gzip.open(filename, mode, encoding="utf-8")
    else:
        return open(filename, mode)


# Ensembl (ENST00000000233.5), FlyBase (FBtr0070000.1) and RefSeq (NM_000005.3) accessions end in a version.
# In other IDs the last number is part of the ID: PacBio (PB.1.1), StringTie (STRG.9.1) or TAIR (AT1G01010.1).
VERSIONED_ID_REGEX = re.compile(r"^(ENS[A-Z]*\d+|FB[a-z]{2}\d+|[A-Z]{2}_\d+)\.\d+$")


def normalize_id(identifier):
    """
    Removes the version from Ensembl, FlyBase and RefSeq accessions (ENSP00000495982.2 -> ENSP00000495982).
    Any other ID is returned unchanged, so PB.1.1 does not become PB.
    """
    match = VERSIONED_ID_REGEX.match(identifier)
    return match.group(1) if match else identifier


TRAILING_NUMBER_REGEX = re.compile(r"\.\d+$")
GTF_PROTEIN_ID_REGEX = re.compile(r'protein_id "([^"]+)"')


def read_gtf_protein_ids(gtf_file):
    """Protein IDs of the CDS lines of a GTF, as the annotation writes them."""
    protein_ids = set()
    with openfile(gtf_file, "rt") as read_handler:
        for line in read_handler:
            if "\tCDS\t" in line:
                match = GTF_PROTEIN_ID_REGEX.search(line)
                if match:
                    protein_ids.add(match.group(1))
    return protein_ids


def annotation_id_forms(fasta_ids, annotation_ids, database):
    """
    Returns {protein FASTA ID: ID as the annotation of its database writes it}, so every protein ID of the
    GFF3 has the form its database publishes: Ensembl GTFs drop the version (ENSP00000000233.5 ->
    ENSP00000000233), RefSeq ones keep it (NP_000005.3) and in TAIR the number is the isoform (AT1G01010.1).
    The form is read from the annotation and not from ID patterns. normalize_id() is for lookups only.

    IDs that are not in the annotation (human Ensembl proteins on scaffolds, which the .chr GTF leaves
    out) follow what the IDs of the same database do: they lose the number only when every ID found in
    the annotation does, and only if no two IDs become the same one.
    """
    same = {fasta_id for fasta_id in fasta_ids if fasta_id in annotation_ids}
    without_number = {fasta_id: TRAILING_NUMBER_REGEX.sub("", fasta_id) for fasta_id in fasta_ids if fasta_id not in same}
    cut = {fasta_id: short_id for fasta_id, short_id in without_number.items() if short_id in annotation_ids}
    unmatched = [fasta_id for fasta_id in without_number if fasta_id not in cut]

    forms = {fasta_id: fasta_id for fasta_id in fasta_ids}
    forms.update(cut)

    n_kept_whole = 0
    if cut and not same:
        times_used = Counter(TRAILING_NUMBER_REGEX.sub("", fasta_id) for fasta_id in fasta_ids)
        for fasta_id in unmatched:
            if times_used[without_number[fasta_id]] == 1:
                forms[fasta_id] = without_number[fasta_id]
            else:
                n_kept_whole += 1

    logging.info(f"{database}: {len(same)} protein IDs are written as in the annotation, {len(cut)} without their "
                 f"final number, as the annotation has them, and {len(unmatched)} are not in the annotation.")
    if unmatched and ((same and cut) or not (same or cut)):
        logging.warning(f"{database}: the annotation and the protein FASTA do not share one form of protein ID "
                        f"({len(same)} IDs are the same, {len(cut)} differ in the final number). The {len(unmatched)} "
                        f"protein IDs that are not in the annotation are written as in the protein FASTA.")
    if n_kept_whole:
        logging.warning(f"{database}: {n_kept_whole} protein IDs that are not in the annotation keep their final "
                        f"number, because without it they are the same ID as another protein.")
    return forms


# A stop codon more than 50 nt before the last exon-exon junction marks the transcript for nonsense-mediated
# decay. Same number as predicted_NMD of SQANTI3.
NMD_MIN_DISTANCE = 50


def bases_to_last_junction(exons, cds_last, strand):
    """
    Number of bases of the spliced transcript between the last base of the CDS and the last exon-exon
    junction. None when the CDS ends in the last exon or the transcript has one exon.

    :param exons: (start, end) of every exon on the genome, both included
    :param cds_last: genomic position of the last CDS base in the sense of the transcript
    :param strand: "+" or "-"
    """
    ordered_exons = sorted(exons, reverse=(strand == "-"))
    distance = None
    for start, end in ordered_exons[:-1]:
        if distance is not None:
            distance += end - start + 1
        elif start <= cds_last <= end:
            distance = end - cds_last if strand == "+" else cds_last - start
    return distance


def predict_nmd(exons, cds_last, strand, stop_in_cds):
    """
    NMD or nonNMD with the 50 nt rule. The distance is measured on the spliced transcript, which is what
    the ribosome reads. The distance on the genome counts the introns between the stop and the last
    junction, so a stop two exons back was always NMD, also when only a few bases of exon separate it
    from the junction (NRAS, NM_002524.5: 43 nt on the transcript, more than 380 on the genome).

    :param stop_in_cds: True when the CDS includes the stop codon (SQANTI3 ORFs). Ensembl and RefSeq CDS
                        lines end before it, so its 3 bases are not part of the distance.
    """
    distance = bases_to_last_junction(exons, cds_last, strand)
    if distance is None:
        return "nonNMD"
    if not stop_in_cds:
        distance -= 3
    return "NMD" if distance > NMD_MIN_DISTANCE else "nonNMD"


class ChromosomeMap:
    """
    Bi-directional, case-insensitive, version-safe chromosome mapper.
    Every name is translated to the naming of target_db ("ensembl" or "refseq"), which must be the naming of
    the genome FASTA the output is compared with. The pipeline always uses the Ensembl genome
    (prepare_ensembl_reference), also in refseq mode. target_names, when given, is the set of names
    that genome has: a translation that is not in it counts as unmapped, so the callers' 50% abort fires
    instead of writing an empty output.
    """
    NCBI_ACCESSION_PREFIXES = ("NC_", "NW_", "NT_", "AC_", "NZ_")

    def __init__(self, mapping_file=None, target_db="ensembl", target_names=None):
        self.target_db = (target_db or "ensembl").lower()
        self.target_names = set(target_names) if target_names is not None else None
        self.forward_map = {}  # Ensembl -> RefSeq
        self.reverse_map = {}  # RefSeq -> Ensembl
        self.norm_map = {}     # normalized_key -> target_canonical
        self.unmapped_log = set()
        self.total_queries = 0
        self.unmapped_queries = 0
        self.mapping_file = mapping_file

        if mapping_file:
            self._load_mapping_file(mapping_file)

    def _normalize_key(self, key):
        if not key:
            return ""
        k = str(key).strip().lower()
        if k.startswith("chr"):
            k = k[3:]
        return k

    def _strip_version(self, key):
        if any(key.startswith(prefix) for prefix in self.NCBI_ACCESSION_PREFIXES):
            return key.split(".")[0]
        return key

    def _load_mapping_file(self, filepath):
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    ensembl, refseq = parts[0].strip(), parts[1].strip()
                    
                    self.forward_map[ensembl] = refseq
                    self.reverse_map[refseq] = ensembl
                    
                    target = refseq if self.target_db == "refseq" else ensembl
                    
                    for val in (ensembl, refseq):
                        self.norm_map[self._normalize_key(val)] = target
                        self.norm_map[self._normalize_key(self._strip_version(val))] = target

    def get(self, query_chr, default=None):
        """Behaves like a dict.get() for seamless drop-in compatibility."""
        if query_chr is None:
            return default

        q = str(query_chr).strip()
        res = self._count(q)
        if res is not None:
            return res

        # Unmapped contig fallback
        return default if default is not None else q

    def _count(self, q):
        """Looks q up and records the result in the failure stats that check_unmapped uses."""
        res = self._lookup(q)
        self.total_queries += 1
        if res is None:
            self.unmapped_queries += 1
            self.unmapped_log.add(q)
        return res

    def _lookup(self, q):
        """Returns the target name of q, or None if it has none."""
        norm_q = self._normalize_key(q)
        # Tier 1: Direct directional lookup
        if self.target_db == "refseq" and q in self.forward_map:
            res = self.forward_map[q]
        elif self.target_db == "ensembl" and q in self.reverse_map:
            res = self.reverse_map[q]
        # Tier 2: Normalized lookup
        elif norm_q in self.norm_map:
            res = self.norm_map[norm_q]
        # Tier 3: Version-stripped normalized lookup
        else:
            res = self.norm_map.get(self._normalize_key(self._strip_version(q)))

        # A name the genome FASTA does not have is as good as unmapped
        if self.target_names is not None and res not in self.target_names:
            return None
        return res

    def __getitem__(self, key):
        # Not counted: callers check "in" first, and that lookup was counted already
        res = self._lookup(str(key).strip())
        if res is None:
            raise KeyError(key)
        return res

    def __contains__(self, key):
        if not key:
            return False
        return self._count(str(key).strip()) is not None

    def check_unmapped(self, caller, max_rate=0.5):
        """
        Stops the run when more than max_rate of the lookups found no name: the chromosome names of the input,
        the mapping table and the genome do not match, and the output would be empty or nearly so.
        Every lookup counts, not every distinct name, so one unmapped chromosome on most lines stops it.
        """
        if self.total_queries and self.unmapped_queries / self.total_queries > max_rate:
            raise RuntimeError(
                f"Critical Error in {caller}: {self.unmapped_queries}/{self.total_queries} chromosome lookups failed "
                f"using mapping file '{self.mapping_file}' (unmapped names: {', '.join(sorted(self.unmapped_log)[:10])}). "
                f"Halting execution to prevent silent data loss."
            )

    def __len__(self):
        return len(self.norm_map)

    def __bool__(self):
        return bool(self.norm_map)

    def __iter__(self):
        return iter(self.keys())

    def values(self):
        """Returns all mapped target chromosome names."""
        return set(self.forward_map.values()) | set(self.reverse_map.values()) | set(self.norm_map.values())

    def keys(self):
        """Returns all recognized source chromosome names."""
        return set(self.norm_map.keys())

    def items(self):
        """Returns (source_chrom, target_chrom) pairs."""
        return [(k, self.get(k)) for k in self.keys()]


def read_chr_ref_acc(filename, target_db="ensembl", target_names=None):
    """
    Reads a file containing Refseq chromosome accessions and Ensembl identifiers.
    Returns a ChromosomeMap that translates names to the naming of target_db.
    """
    if not filename:
        return ChromosomeMap(target_db=target_db, target_names=target_names)
    return ChromosomeMap(mapping_file=filename, target_db=target_db, target_names=target_names)


def query_biomart_with_retry(dataset, attributes, layer_name="biomart", max_retries=5, initial_delay=5, disable_hint=True):
    """
    Queries pybiomart Dataset with up to max_retries attempts using exponential backoff.
    If all attempts fail, raises RuntimeError with explicit instructions on how to disable the layer,
    or a saved table can be used instead.

    :param disable_hint: add the advice on turning the layer off; False when the caller is not a layer.
    """
    import time
    import logging

    for attempt in range(1, max_retries + 1):
        try:
            logging.info(f"[{layer_name}] Querying BioMart (Attempt {attempt}/{max_retries})...")
            result = dataset.query(attributes=attributes)
            # When Ensembl is down, pybiomart reads its "Service unavailable" HTML page as a one-column
            # table instead of raising. Check the columns so that page is retried like any error.
            if result.shape[1] != len(attributes):
                raise ValueError(f"BioMart answered with {result.shape[1]} column(s) instead of {len(attributes)}, "
                                 f"probably an error page: {str(result.columns[0])[:80]}")
            return result
        except Exception as err:
            if attempt == max_retries:
                # The Reactome layer is turned off with the reactome key; layer_reactome=no has no effect
                config_key = "reactome" if layer_name == "layer_reactome" else layer_name
                message = f"Fatal Error in '{layer_name}': BioMart query failed after {max_retries} attempts due to server/network error ({err})."
                if disable_hint:
                    message += (f" To run IsoAnnot without BioMart, deactivate this layer by passing '--config {config_key}=no' to isoannot.sh, "
                                f"or set '{config_key}: no' in your species config.yaml.")
                raise RuntimeError(message)
            delay = initial_delay * (2 ** (attempt - 1))
            logging.warning(f"[{layer_name}] BioMart attempt {attempt}/{max_retries} failed ({err}). Retrying in {delay}s...")
            time.sleep(delay)


# BioMart attributes saved by save_biomart_tables.py, and the column names BioMart gives them.
# layer_go uses the GO columns and layer_reactome the gene ID and name pairs.
BIOMART_TABLE_ATTRIBUTES = ["ensembl_gene_id", "external_gene_name", "go_id", "name_1006", "namespace_1003"]
BIOMART_TABLE_COLUMNS = ["Gene stable ID", "Gene name", "GO term accession", "GO term name", "GO domain"]


def read_biomart_table(filename, layer_name="biomart"):
    """
    Reads a table written by save_biomart_tables.py: the one the get_biomart_table rule downloads, or the
    one of the biomart_table config key. Empty cells become NA, as in a BioMart answer.
    """
    import logging
    import pandas as pd

    table = pd.read_csv(filename, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    missing = [col for col in BIOMART_TABLE_COLUMNS if col not in table.columns]
    if missing:
        raise ValueError(f"{filename} is not a saved BioMart table: columns {', '.join(missing)} are missing. "
                         f"Write it with scripts/save_biomart_tables.py.")
    logging.info(f"[{layer_name}] Using the saved BioMart table {filename} ({len(table)} rows)")
    return table


def merge_fasta_dicts(fasta_files):
    """
    Merges several fasta containing dictionaries.

    Args: 
        fasta_files (dict): dictionary of dictionaries containing fasta files
            {tag of fasta_file:{fasta_key:fasta:_value, fasta_key:fasta_value}}

    Returns: merged_fasta (dict): dictionary containing merged fasta files
    """
    merged_fasta = defaultdict(lambda: defaultdict(list))
    for tag, fasta_file in fasta_files.items():
        for fasta_key, fasta_value in fasta_file.items():
            merged_fasta[fasta_key][tag].extend(fasta_value)
    return merged_fasta


def natural_key(text):
    """
    Sort key that compares the numbers inside an ID as numbers: Q9HBE1-2 < Q9HBE1-10, PB.2.1 < PB.10.1.
    Used to pick protein IDs in a fixed order (GitHub #10).
    """
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", text)]


def get_consecutive_parts(position_list):
    """
    From a list of positions, it returns a list of consecutive positions
    (eg. CDS consecutive positions). Each element of the list contains a list 
    with the start and end positions of the consecutive part.
    Args: 
        position_list (list): list of positions.

    Returns: 
        loc_range (list): 
    
    """
    # TODO: sort always?
    for k, g in itertools.groupby(enumerate(sorted(position_list)), lambda t: t[0] - t[1]):
        loc_range = list(map(itemgetter(1), g))
        yield loc_range


def get_summary_parts(position_list):
    """
    Generates a dictionary with the start, end and length of each 
    consecutive part.
    Args: 
        position_list (list): list of positions.

    Returns: 
        output dict (dict): containing start, end and length of each
        consecutive part.
    """
    output_dict = {
        "start": [],
        "end": [],
        "len": 0
    }

    for loc_range in get_consecutive_parts(position_list):
        output_dict["start"].append(loc_range[0])
        output_dict["end"].append(loc_range[-1])
        output_dict["len"] += loc_range[-1] - loc_range[0] + 1

    return output_dict


# Function is no longer needed, we are using refseq and ensembl GTF (only needed when using GFF)
# def _remove_stop_codon(df_group):  
#     """
    
#     """
#     # For positive strands remove 3 positions from latest end position
#     if "+" in df_group["Strand"].values:
#         df_group.loc[df_group["End"].idxmax(), "End"] -= 3
#     # For negative strands increase 3 positions from the first start position
#     elif "-" in df_group["Strand"].values:
#         df_group.loc[df_group["Start"].idxmin(), "Start"] += 3

#     return df_group

# def read_feature_from_gtf(gtf_file:str, feature:str, attribute:str, chr_ref=None, match_prefix=""):
#     """
#     Gets info from GTF file

#     Args: 
#         gtf_file (str): name of gtf file
#         feature (str): feature you want to retrieve
#         attribute (str):
#         chr_ref ():
#         match_prefix(str):
#         correct_stop_codon (boolean):
    
#     Returns:
#         gtf_contents (dict): dictionary with gtf content
#     """
#     Returns pyranges as dataframe
#     gtf_contents = pyranges.read_gtf(gtf_file, as_df=True, duplicate_attr=True)
#     Keep 1-based as original GTF
#     gtf_contents['Start'] = gtf_contents['Start'].apply(lambda x: x + 1)

#     Keep gtf_contents only if it contains the feature we are searching for and has attributes
#     gtf_contents = gtf_contents[(gtf_contents["Feature"] == feature) & (gtf_contents[attribute].notna())]

#     if chr_ref:
#         dict_table = pd.read_csv(chr_ref, sep="\t", comment="#", header=None)
#         Append prefix
#         dict_table.iloc[:,0] = match_prefix + dict_table.iloc[:,0].astype(str)
#         replace_dict = dict(zip(dict_table.iloc[:,1], dict_table.iloc[:,0]))
#         gtf_contents["Chromosome"] = gtf_contents["Chromosome"].map(replace_dict) 

#     return gtf_contents
