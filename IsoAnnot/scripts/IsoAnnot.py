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


def read_chr_ref_acc(filename):
    """
    Reads a file containing Refseq chromosome accessions
    and their corresponding Ensembl identifyer.

    Args:
        filename (str): name of file.
    Returns:
        chr_accession_dict (dict): {refseq_acc:ensembl_acc}
    
    """
    refseqChrom = {}
    with open(filename, 'r') as infile:
        for line in infile:
            if line[0] != "#":
                refseqChrom[line.split("\t")[1].strip()] = line.split("\t")[0].strip()

    return refseqChrom


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
