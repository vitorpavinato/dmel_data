"""
This modeule contains a set of functions to manipulate SNP data stored 
in a pandas DataFrame. 
"""

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

    # Raise error for empty tst_df and bed_df
    if tsv_df.empty:
        raise ValueError("tsv_df is empty")

    if bed_df.empty:
        raise ValueError("bed_df is empty")

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


# Impute missing haplotypes
def impute_missing_haplotypes(input_df: DataFrame, max_number_haplotypes: int) -> DataFrame:
    """
    This function imputes missing haplotypes by replacing the missing values, 
    that is the difference between the total count and the maximum total count.
    This function replaces complete the number of genotypes only for the common
    allele.
    """

    # Raise error for an empty DataFrame
    if input_df.empty:
        raise ValueError("DataFrame is empty")

    # Copy the input dataframe
    df = input_df.copy()

    for index, row in df.iterrows():  # Iterate over the rows of the DataFrame
        if row['totalcount'] < max_number_haplotypes:
            allele_counts = (row['refcount'], row['altcount'])  # Access counts for the current row
            idx_common_allele = allele_counts.index(max(allele_counts))
            if idx_common_allele == 0:
                df.at[index, 'refcount'] += max_number_haplotypes - row['totalcount']  # Update refcount for the current row
            else:
                df.at[index, 'altcount'] += max_number_haplotypes - row['totalcount']  # Update altcount for the current row

    return df


# Convert a given pandas DataFrame to a dictionary of SNPs counts
def create_snp_total_counts_dict(df: DataFrame, ) -> dict[int, list[int]]:
    """
    Covert a given pandas DataFrame containing SNPs to a
    dictionary of SNP total counts. The key of the dictionary
    is the each unique SNP total counts. Each value is a
    list of altcounts for each total counts. It expect a df with
    "refcount" and "altcount" fields. It returns a sorted dictionary
    of lists of altcounts for each total counts.
    """

    # Raise error for an empty DataFrame
    if df.empty:
        raise ValueError("DataFrame is empty")

    # Initialize the snps dictionary
    snp_total_counts_dict = {}

    # Fill up the SNPs dictionary using the table data
    for _, row in df.iterrows():
        ref_counts = row["refcount"]
        alt_counts = row["altcount"]
        total_counts = ref_counts + alt_counts

        if total_counts not in snp_total_counts_dict:
            snp_total_counts_dict[total_counts] = []
        snp_total_counts_dict[total_counts].append(alt_counts)

    snp_total_counts_dict_sorted = dict(sorted(snp_total_counts_dict.items()))

    return snp_total_counts_dict_sorted
