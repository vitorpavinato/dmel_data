"""
This module contains a set of functions to create 
site-frequency spectrums (SFSs) and manipulate SFSs.
"""

from scipy.stats import hypergeom


# Get the a SFS from each of SNP total counts in dict
def create_sfs_dict(
    snp_counts_dict: dict[int, list[int]]
) -> dict[int, list[int]]:
    """
    Generate a SFS from each sample sizes (total counts) in dict.
    The SFS is a list of counts for each bin in the SFS. It expects
    a sorted dict of lists of altcounts for each total counts.
    It returns a sorted dict of SFS (which is a list of bin counts) for each
    total counts.
    """

    # Raise error for an empty dictionary
    if not snp_counts_dict:
        raise ValueError("Dictionary is empty")

    # Sort the dictionary of lists of altcounts for each total counts
    snp_counts_dict = dict(sorted(snp_counts_dict.items()))

    # Create an empty dict to save each popsize SFS
    # list of values for each key popsize
    snp_total_counts_dict = {}

    for key, value in snp_counts_dict.items():
        snp_total_counts_dict[key] = [0] * (key + 1)
        for alt_counts in value:
            snp_total_counts_dict[key][alt_counts] += 1

    return snp_total_counts_dict


# Downsampling SFSs
def downsample_sfs(
    original_sfs: list[int],
    original_size: int,
    sample_size: int
) -> list[int | float]:
    """
    Project a distribution an unfolded or folded site-frequency spectrum  
    to a sample size m < n.
    :param original_sfs: List of counts for each bin in the original SFS.
    Expects an SFS as a list of length sample size (n) + 1,
    with position i for i gene copies;
    :param popsize: Total number of gene copies in the original distribution.
    :param sampsize: Total number of gene copies in the
    downsampled distribution (m < n).
    :return: List of expected counts for each bin in the
    projected distribution.
    """

    # Raise error for an empty SFS
    if not original_sfs:
        raise ValueError("SFS is empty")

    # Create an empty sample SFS with sample_size + 1 bins
    sample_sfs = [0]*(sample_size+1)

    for pi, count in enumerate(original_sfs):
        for si in range(sample_size+1):
            prob = hypergeom.pmf(si, original_size, pi, sample_size)
            sample_sfs[si] += prob * count

    return sample_sfs


# Fold an unfolded SFS
def fold_sfs(
    sfs: list[int | float]
) -> list[int | float]:
    """
    Fold an unfolded SFS
    """

    # Raise error for an empty SFS
    if not sfs:
        raise ValueError("SFS is empty")

    # Get the middle index
    mid_idx = int(len(sfs)/2)

    if len(sfs) % 2 > 0:
        
        mid_value = sfs.pop(mid_idx)
    else:
        mid_value = None

    # Devide the SFS in two parts:
    # First part before the mid-index, the second part
    # after the mid-index. The two list should no include
    # the mid-index. Then revert the second list.
    before_mididx = sfs[0:mid_idx]
    after_mididx = sfs[mid_idx:]

    rev_after_mididx = after_mididx.copy()
    rev_after_mididx.reverse()

    fsfs = [x + y for x, y in zip(before_mididx, rev_after_mididx)]

    if mid_value is not None:
        fsfs.append(mid_value)

    return fsfs