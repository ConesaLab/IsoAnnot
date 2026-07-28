"""
IsoAnnot functions
Modified by Alessandra Martinez
"""
import gzip, itertools
from collections import defaultdict
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


class ChromosomeMap:
    """
    Bi-directional, case-insensitive, version-safe chromosome mapper.
    Supports Ensembl-leading (default), RefSeq-leading, or identity mapping mode.
    """
    NCBI_ACCESSION_PREFIXES = ("NC_", "NW_", "NT_", "AC_", "NZ_")

    def __init__(self, mapping_file=None, leading_db="ensembl"):
        self.leading_db = (leading_db or "ensembl").lower()
        self.forward_map = {}  # Ensembl -> RefSeq
        self.reverse_map = {}  # RefSeq -> Ensembl
        self.norm_map = {}     # normalized_key -> target_canonical
        self.unmapped_log = set()
        self.total_queries = 0

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
                    
                    target = refseq if self.leading_db == "refseq" else ensembl
                    
                    for val in (ensembl, refseq):
                        self.norm_map[self._normalize_key(val)] = target
                        self.norm_map[self._normalize_key(self._strip_version(val))] = target

    def get(self, query_chr, default=None):
        """Behaves like a dict.get() for seamless drop-in compatibility."""
        if query_chr is None:
            return default

        q = str(query_chr).strip()
        self.total_queries += 1

        # Tier 1: Direct directional lookup
        if self.leading_db == "refseq" and q in self.forward_map:
            return self.forward_map[q]
        if self.leading_db == "ensembl" and q in self.reverse_map:
            return self.reverse_map[q]

        # Tier 2: Normalized lookup
        norm_q = self._normalize_key(q)
        if norm_q in self.norm_map:
            return self.norm_map[norm_q]

        # Tier 3: Version-stripped normalized lookup
        stripped_q = self._normalize_key(self._strip_version(q))
        if stripped_q in self.norm_map:
            return self.norm_map[stripped_q]

        # Unmapped contig fallback
        self.unmapped_log.add(q)
        return default if default is not None else q

    def __getitem__(self, key):
        res = self.get(key, default=None)
        if res is None and key not in self.norm_map:
            raise KeyError(key)
        return res

    def __contains__(self, key):
        if not key:
            return False
        q = str(key).strip()
        return (q in self.forward_map or q in self.reverse_map or
                self._normalize_key(q) in self.norm_map or
                self._normalize_key(self._strip_version(q)) in self.norm_map)


def read_chr_ref_acc(filename, leading_db="ensembl"):
    """
    Reads a file containing Refseq chromosome accessions and Ensembl identifiers.
    Returns a ChromosomeMap instance.
    """
    if not filename:
        return ChromosomeMap(leading_db=leading_db)
    return ChromosomeMap(mapping_file=filename, leading_db=leading_db)


def query_biomart_with_retry(dataset, attributes, layer_name="biomart", max_retries=5, initial_delay=5):
    """
    Queries pybiomart Dataset with up to max_retries attempts using exponential backoff.
    If all attempts fail, raises RuntimeError with explicit instructions on how to disable the layer.
    """
    import time
    import logging

    for attempt in range(1, max_retries + 1):
        try:
            logging.info(f"[{layer_name}] Querying BioMart (Attempt {attempt}/{max_retries})...")
            return dataset.query(attributes=attributes)
        except Exception as err:
            if attempt == max_retries:
                raise RuntimeError(
                    f"Fatal Error in '{layer_name}': BioMart query failed after {max_retries} attempts due to server/network error ({err}). "
                    f"To run IsoAnnot without BioMart, deactivate this layer by passing '--config {layer_name}=no' to isoannot.sh, "
                    f"or set '{layer_name}: no' in your species config.yaml."
                )
            delay = initial_delay * (2 ** (attempt - 1))
            logging.warning(f"[{layer_name}] BioMart attempt {attempt}/{max_retries} failed ({err}). Retrying in {delay}s...")
            time.sleep(delay)


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
