'''
A set of function to create context dependent site-frequency spectrums (SFSs).
SNPs in the SFS might be constrained by distance, by similar mutational context
(here defined as the target SNP mutation plus and minos 1 nucleoted: AA/TC for
example). The functions included here help identify the distance contraint,
by calculating all the distances among synonymous and non-synonymous SNPs for
example ("calculate_distances"). It allows also identify blocks of consecutive
SNPs that might createambiguous mutational context
("find_consecutive_positions"). The remaining functions are part of the
mutational context pairing of the original pipeline, kept here for the record.
'''
from typing import Collection
import numpy as np
from pandas import DataFrame
from .sfs_utils import downsample_sfs, fold_sfs


def calculate_distances(
    df: DataFrame, type1: str, type2: str
) -> Collection[int | float]:
    """
    Giving a dataFrame, and two mutation types,
    this function returns the absolute distance
    between the type1 and type2. The dataFrame
    should have "effect" and "pos" fields.
    """

    # Filter rows for the specified types
    type1_df = df[df["effect"] == type1]
    type2_df = df[df["effect"] == type2]

    # Initialize a list to store distances
    distances = []

    # Iterate through each entry of type1
    for _, row1 in type1_df.iterrows():
        pos1 = row1["pos"]

        # Iterate through each entry of type2
        for _, row2 in type2_df.iterrows():
            pos2 = row2["pos"]

            # Calculate the absolute difference in positions
            distance = abs(pos1 - pos2)

            # Append the distance to the list
            distances.append(distance)

    return distances


def find_consecutive_positions(positions: list[int]) -> list[bool]:
    """
    This function find blocks of consecutive SNPs and
    returns False if a SNPs is part of a block of SNPs
    and True otherwise. This is useful to find and remove
    positions that might make ambiguous mutational
    contexts. Because of the nature of each of mutational
    context were defined based on the reference and
    alternative alleles.
    """
    pos = positions

    # Initialize a list to store consecutive blocks
    consecutive_blocks = []

    # Initialize variables for the current block
    current_block = [pos[0]]

    # Iterate through the positions starting from the second element
    for i in range(1, len(pos)):
        # Check if the current position is consecutive to the previous one
        if pos[i] - pos[i - 1] == 1:
            current_block.append(pos[i])
        else:
            # If not consecutive, start a new block
            consecutive_blocks.append(current_block)
            current_block = [pos[i]]

    # Append the last block
    consecutive_blocks.append(current_block)

    positions_to_keep = []
    for block in consecutive_blocks:
        if len(block) == 1:
            positions_to_keep.append(True)
        else:
            for i in block:
                positions_to_keep.append(False)

    return positions_to_keep


def set_snp_sequence_category(df: DataFrame) -> list[str]:
    """
    This define the mutational context each SNP in a dataFrame.
    Create a column containing the sequence category for each SNP
    in a DataFrame. The df should have, "ref", "alt" and "altflank"
    fields.
    """
    list_refalleles = df['ref']
    list_altalleles = df['alt']
    list_altflanks = df['altflank']
    
    sequence_categories = [f"{altflank[3]}{ref}/{alt}{altflank[5]}" for ref, alt, altflank  in zip(list_refalleles, list_altalleles, list_altflanks)]
    return sequence_categories


def create_sequence_categories_dic() -> dict[int, str]:
    """
    This creates a look-up dictionary for the sequence
    category. The dictionary then needs to be provided to
    "create_snp_dict" or "create_snp_dict_wrapper" in order
    to make a SNP dataFrame a dictionary of SNP organized by
    their SNP category.
    """

    # Possible SNPs:
    snp_types = [["A", "C"],
                 ["A", "G"],
                 ["A", "T"],
                 ["C", "G"],
                 ["C", "T"],
                 ["G", "T"]]

    # List of DNA nucleotides
    nucleotides = ["A", "C", "G", "T"]

    # Create an empty dictionary to save
    # the SNPs categories
    sequence_categories_dict = {}

    # Start the first dictionary key
    key_number = 1

    for snp_pairs in snp_types:
        for j in nucleotides:
            for i in nucleotides:
                key = key_number
                value = [
                         i + snp_pairs[0] + "/" + snp_pairs[1] + j,
                         i + snp_pairs[1] + "/" + snp_pairs[0] + j
                         ]
                sequence_categories_dict[key] = value

                key_number += 1
                # # Set the reciprocal
                # key = key_number+1
                # value = i + snp_pairs[1] + "/" + snp_pairs[0] + j
                # sequence_categories_dict[key] = value
                # key_number += 2

    return sequence_categories_dict


def create_snp_dict(
    df: DataFrame,
    sequence_categories_dict: dict[int, str]
) -> dict[int, dict[int, list[int, int]]]:
    """
    Covert a given pandas DataFrame containing SNPs
    to a dictionary. It needs a sequence category dictionary
    created with "create_sequence_categories_dic". The df
    should have "refcount", "altcount" and "sequence_category"
    fields.
    """

    # Initialize the snps dictionary
    snps_dict = {}

    # Fill up the SNPs dictionary using the table data
    for _, row in df.iterrows():
        pos = row['pos']
        ref_counts = row['refcount']
        alt_counts = row['altcount']
        total_counts = ref_counts + alt_counts
        sequence_category = row['sequence_category']

        # Find the corresponding numeric key(s) from the lookup dictionary
        numeric_keys = [
            key for key, value in sequence_categories_dict.items()
            if sequence_category in value
            ]

        # Add the data to the selected_snps dictionary for each numeric key
        for numeric_key in numeric_keys:
            # numeric_key_int = int(numeric_key)
            if numeric_key not in snps_dict:
                snps_dict[numeric_key] = {}
            snps_dict[numeric_key][pos] = [alt_counts, total_counts]

    return snps_dict


def create_snp_dict_wrapper(
    df: DataFrame,
    sequence_categories_dict: dict[int, str]
) -> dict[int, dict[int, list[int, int]]]:
    """
    Wrapper for create_snp_dict. This wrapper allows to process different
    effect mutations: introns, non-synonymous and synonymous.
    It expect the "effect" field to be one of "INTRON",
    "NON_SYNONYMOUS_CODING".
    """
    introns_df = df[df['effect'] == "INTRON"]
    nonsyns_df = df[df['effect'] == "NON_SYNONYMOUS_CODING"]
    syns_df = df[df['effect'] == "SYNONYMOUS_CODING"]

    introns_dict = create_snp_dict(introns_df, sequence_categories_dict)
    nonsyns_dict = create_snp_dict(nonsyns_df, sequence_categories_dict)
    syns_dict = create_snp_dict(syns_df, sequence_categories_dict)

    return introns_dict, nonsyns_dict, syns_dict


def find_closest(lst: list[int], k: int) -> tuple[int, int] | None:
    """
    Find the closest number to k in a list of numbers.
    This functions is part of the find_close_snp_pairs.
    """
    # lst.sort()
    closest_num_idx = 0
    closest_num = lst[0]
    if len(lst) > 1:
        for idx, num in enumerate(lst):
            if abs(num - k) < abs(closest_num - k):
                closest_num = num
                closest_num_idx = idx
            if num > k:
                break
    return closest_num, closest_num_idx


def simple_find_closest(lst: list[int], k: int) -> tuple[int, int] | None:
    """
    This simply take the first item in one list as the closest
    item to the other list. The list of position is sorted,
    so the first item should be the closest item to the first
    other list.
    """
    # lst.sort()
    closest_num_idx = 0
    closest_num = lst[0]
    return closest_num, closest_num_idx


def find_closest_snp_pairs(
    dict1: dict[int, dict[int, list[int, int]]],
    dict2: dict[int, dict[int, list[int, int]]]
) -> dict[int, dict[int, list[int, int]]]:
    """Find pairs of closest SNPs from two distinct
    dictionaries. The first dictionary must be the one
    with less SNPs (usually introns or any other neutral
    control)
    """

    # Check if if both dictionaries are not empty
    if dict1 and dict2:
        print("Both dictionaries are not empty!")
    else:
        raise ValueError("One or both dictionaries are empty!")

    # Two empty dictionaries to save each SNP
    # in a pair found in both input dictionaries
    paired_snps_dict1 = {}
    paired_snps_dict2 = {}

    # Transverse each sequence category (seq_key) in the first
    # dictionary, finding the same seq_key in the second dictionary
    for seq_key in dict1:
        if seq_key in dict2:
            dict1_positions = dict1[seq_key]
            dict2_positions = dict2[seq_key]
            list_dict2_positions = list(dict2_positions)

            # When there is a seq_key in both, transvers the list of 
            # positions in the first, finding a closes SNPs in the second.
            # Remove the kept SNPs from the second dictionary list to 
            # avoid repeting the same SNP. Leave the loop when the list of
            # SNPs in the second dictionary seq_key is exhausted.

            for pos in dict1_positions:
                # closest_pos, closest_pos_idx = find_closest(list_dict2_positions, pos)
                closest_pos, closest_pos_idx = simple_find_closest(list_dict2_positions, pos)

                paired_snps_dict1.setdefault(seq_key, {})[pos] = dict1_positions[pos]
                paired_snps_dict2.setdefault(seq_key, {})[closest_pos] = dict2_positions[closest_pos]

                # Remove the selected position from the list
                # This works like a FIFO algorithm
                list_dict2_positions.pop(closest_pos_idx)
                if len(list_dict2_positions) == 0:
                    break

    return paired_snps_dict1, paired_snps_dict2


# Maybe remove this function
def sfs_from_snp_dict(
    snp_dict: dict[int, dict[int, list[int, int]]],
    min_sample_size: int,
    max_sample_size: int,
    folded: bool,
) -> list[int | float]:
    """
    Get the SFS from each sample size in dict. The input
    dictionary should have the following hierarchy:
    sequence_category: position: [alt_count, total_count].
    """

    # Create an empty dict to save each popsize SFS
    # list of values for each key popsize
    sample_size_dict = {}

    for i in range(min_sample_size, max_sample_size + 1):
        sample_size_dict[i] = [0] * (i+1)

    # Sequence category -wise loop to fill the
    # dictionary of popsize SFS
    for seq_key in snp_dict:
        dict_pos = snp_dict[seq_key]
        for pos in dict_pos:
            snp_counts = dict_pos[pos]
            sample_size_dict[snp_counts[1]][snp_counts[0]] += 1

    # Convert the sample_size_dict to a list and apply
    # downsampling in the SFS for higher values than
    # min_sample_size
    sfs_list = []
    for sample_popsize, sfs in sample_size_dict.items():
        if sample_popsize == min_sample_size:
            sfs_list.append(sfs)
        else:
            ds = downsample_sfs(sfs, sample_popsize, min_sample_size)
            sfs_list.append(ds)

    # Conver the list of SFS to an np.array
    sfs_array = np.array(sfs_list)

    # Sum column-wise to get the final SFS
    sfs = list(np.sum(sfs_array, 0))

    if folded:
        return fold_sfs(sfs)
    return sfs
