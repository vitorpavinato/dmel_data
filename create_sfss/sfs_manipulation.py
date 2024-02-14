"""This module contains a collection of functions
used in the jupyter notebook. It complies to modular coding.
"""

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from pandas import DataFrame


# Create a .BED file with only short introns
def filter_short_introns_from_bed(
    input_bed: str, output_bed: str, chrom_list: list[str],
    short_intron_size: int = 86, trailling_size: int = 8
) -> DataFrame:
    """
    Create a .BED files with only short introns intervals
    """

    short_introns = []

    with open(input_bed, "r", encoding="utf-8") as inbed, open(output_bed, "w", encoding="utf-8") as outbed:
        for linen, line in enumerate(inbed):
            line = line.rstrip("\n")
            if (line.startswith("##") or line.startswith("#") or line.startswith(">") or
                line.startswith("A") or line.startswith("T") or line.startswith("C") or line.startswith("G")):
                continue

            fields = [splits for splits in line.split("\t") if splits != ""]

            if len(fields) <= 1:
                continue

            # Unpack FIELDS items
            chrom, estart, eend, *_ = fields

            # Cast estarts and eends to integer
            estart = int(estart)
            eend = int(eend)

            # Process only chroms in the list
            if any(x == chrom for x in chrom_list):
                if eend - estart < short_intron_size:
                    estart_trimmed = estart + trailling_size
                    eend_trimmed = eend - trailling_size
                    short_intron = [chrom, (estart_trimmed-1), eend_trimmed]
                    short_intron.extend(_)
                    outbed.write('\t'.join(str(item) for item in short_intron) + "\n")
                    short_introns.append(short_intron)

    return DataFrame(short_introns)


# Filter SNPs by interval
def filter_snps_by_interval(
    tsv_df: DataFrame, bed_df: DataFrame
) -> DataFrame:
    """
    Giving a dataFrame containing SNPs, and another,
    dataFrame containing intervals, this function returns
    another dataFrame containing only SNPs that fall
    within the intervals.
    :param tsv_df assumes that chrom name (same as in bed_df)
    is in column with index 0.
    :param bed_df assums a bed-like file with at least 3 columns:
    1. with chrom name (same as in tsv_df);
    2. with starting coordinates for the interval;
    3. with ending coordinates for the interval.
    """

    # Get unique chromosomes from the SNP DataFrame
    unique_chromosomes = tsv_df.iloc[:, 0].unique()

    # Filter bed DataFrame for intervals corresponding to SNP chromosomes
    relevant_intervals = bed_df[bed_df.iloc[:, 0].isin(unique_chromosomes)]

    # Get the starts and ends from the bed DataFrame
    # Get the starts and ends from the relevant intervals
    bed_starts = relevant_intervals.iloc[:, 1].tolist()
    bed_ends = relevant_intervals.iloc[:, 2].tolist()

    # Check if each SNP position falls within any interval
    filtered_df = tsv_df[
        tsv_df.iloc[:, 1].apply(
            lambda pos: any(
                start <= pos <= end for start, end in zip(bed_starts, bed_ends)
            )
        )
    ]
    return filtered_df


# Convert a given pandas DataFrame to a dictionary of SNPs counts
def snp_totalcounts_dict(
    df: DataFrame,
) -> dict[int, list[int, int]]:
    """
    Covert a given pandas DataFrame
    containing SNPs to a dictionary based on
    the SNP total counts. It expect a df with
    "refcount" and "altcount" fields.
    """

    # Initialize the snps dictionary
    snps_counts_dict = {}

    # Fill up the SNPs dictionary using the table data
    for _, row in df.iterrows():
        ref_counts = row["refcount"]
        alt_counts = row["altcount"]
        total_counts = ref_counts + alt_counts

        if total_counts not in snps_counts_dict:
            snps_counts_dict[total_counts] = []
        snps_counts_dict[total_counts].append(alt_counts)

    return snps_counts_dict

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
    sample = [0]*(sample_size+1)

    for pi, count in enumerate(original_sfs):
        for si in range(sample_size+1):
            prob = hypergeom.pmf(si, original_size, pi, sample_size)
            sample[si] += prob * count

    return sample


# Fold an unfolded SFS
def fold_sfs(
    sfs: list[int | float]
) -> list[int | float]:
    """Fold an unfolded SFS"""

    # Find the indexs that divides the sfs in two:
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


# Get the SFS from each sample sizes (total counts) in dict
# Get a combined SFS from a list of
def sfs_from_counts_dict(
    snps_counts_dict: dict[int, list[int]],
    min_size: int,
    max_size: int,
    folded: bool,
) -> list[int | float]:
    """
    Get the SFS from each sample sizes (total counts) in dict.
    Then combine the SFSs.
    """

    # Create an empty dict to save each popsize SFS
    # list of values for each key popsize
    sizes_dict = {}

    for i in range(min_size, max_size + 1):
        sizes_dict[i] = [0] * (i + 1)

    # Sequence category -wise loop to fill the
    # dictionary of popsize SFS
    for si, list_alt_counts in snps_counts_dict.items():
        for alt_counts in list_alt_counts:
            sizes_dict[si][alt_counts] += 1

    # Convert the sizes_dict to a list and apply
    # downsampling in the SFS for higher values than
    # min_size
    sfs_list = []
    for si, sfs in sizes_dict.items():
        if si == min_size:
            sfs_list.append(sfs)
        else:
            ds_sfs = downsample_sfs(sfs, si, min_size)
            sfs_list.append(ds_sfs)

    # Conver the list of SFS to an np.array
    sfs_array = np.array(sfs_list)

    # Sum column-wise to get the final SFS
    sfs = list(np.sum(sfs_array, 0))

    if folded:
        return fold_sfs(sfs)
    return sfs


def get_raw_sfs_from_counts_dict(
    snps_counts_dict: dict[int, list[int]],
    min_size: int,
    max_size: int,
    folded: bool,
) -> dict[int, list[int]]:

    # Create an empty dict to save each popsize SFS
    # list of values for each key popsize
    sizes_dict_sfs = {}

    for i in range(min_size, max_size + 1):
        sizes_dict_sfs[i] = [0] * (i + 1)

    # Sequence category -wise loop to fill the
    # dictionary of popsize SFS
    for si, list_alt_counts in snps_counts_dict.items():
        for alt_counts in list_alt_counts:
            sizes_dict_sfs[si][alt_counts] += 1

    return sizes_dict_sfs


def main() -> None:
    """Test code"""

    # Example usage for filter_snps
    # Generate random data for vcf DataFrame
    np.random.seed(42)  # Setting seed for reproducibility
    tsv = pd.DataFrame(
        {
            "chrom": "chr4",
            "pos": np.sort(np.random.randint(1, 10001, size=20)),
            "effect": np.random.choice(["intron", "non-synonymous"], size=20),
        }
    )

    # Sort the vcf file
    tsv = tsv.sort_values(by=["pos"])

    # Generate random data for bed DataFrame
    bed = pd.DataFrame(
        {
            "chrom": "chr4",
            "starts": [1000, 4000, 6000, 7000, 9000],
            "ends": [1500, 4300, 6999, 7100, 9500],
            "interval": ["intron1", "intron2", "intron3", "intron4", "intron5"],
        }
    )

    # Display the DataFrames
    print("vcf DataFrame:")
    print(tsv)
    print("\nbed DataFrame:")
    print(bed)

    # Filter vcf SNPs
    result_filtered = filter_snps_by_interval(tsv, bed)

    # Display the resulting filtered DataFrame
    print(result_filtered)

    # Example usage for project_sfs_v*
    # A small SFS
    #sfs = [0, 10, 9, 4, 5, 6, 3, 2]
    sfs = [0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    print(sfs, sum(sfs))

    n = 11  # Total number of gene copies
    m = 10  # Sample size m < n

    downsampled_sfs = downsample_sfs(sfs, n, m)
    print(downsampled_sfs, sum(downsampled_sfs))

    # Test the sfs_from_counts_dict function
    test_count_dict = {
                        11: [1, 2, 3, 4, 5],
                        12: [1, 2, 3, 4, 5]

    }
    result = sfs_from_counts_dict(test_count_dict, min_size=11, max_size=12, folded=False)
    print(result, sum(result))


if __name__ == "__main__":
    main()
