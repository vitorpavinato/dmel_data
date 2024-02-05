"""
Sript to simplify SNPEff annotation
For the moment it only takes the first term
It might be redundant as SNPEff has this option
But the idea is to have a more probabilistic approach
in selecting annotated terms.
"""
import os
import argparse
import sys


# FUNCTIONS
def simplify_snpeff(file: str, outfile: str,
                    keeponlyterms: bool = False, method: str = "First") -> None:

    """
    This function simplifies SNPEff entries of a VCF file.
    For the moment, method don't do anything.
    """

    path, filename = os.path.split(file)

    # Define the basename for the outputs from the filename
    basename, _ = filename.strip().split('.vcf')

    # Define input and output files
    inputfile = file

    if outfile is None:
        outputfile = path + basename + "_simplified.vcf"
    else:
        outputfile = outfile

    # Copy the header from the inputfile to the ouputfile
    with open(inputfile, "r", encoding="utf-8") as input_file, open(outputfile, "w", encoding="utf-8") as output_file:
        for line in input_file:
            if (line.startswith("##") or line.startswith("#")):
                output_file.write(line)

    # Define annotation effects terms we care about
    relevant_effect_terms = ['INTRON', 'SYNONYMOUS_CODING', 'NON_SYNONYMOUS_CODING', 'INTERGENIC']

    # Open the input file in read mode and output file in write mode
    with open(inputfile, "r") as input_file, open(outputfile, "a") as output_file:
        # For each line, breakdow the SNPEff info field to retain only the effect entries
        # Iterate through the lines in the input file
        for line in input_file:
            if (line.startswith("##") or line.startswith("#")):
                continue
        
            # Split the line into fields using tab as the delimiter
            fields = line.strip().split('\t')

            # Get the ANN field (assuming it's always the eighth field)
            info_field = fields[7]

            # Split the ANN field by commas to separete the SNPEff annotation
            altcount, refcount, snpeffs, *lo = info_field.split(';')

            if (snpeffs == "ReverseComplementedAlleles"):  # check later if this is still necessary
                _, snpeffann, = lo[0].split("EFF=")
            else:
                _, snpeffann, = snpeffs.split("EFF=")

            # Separete each SNPEff entry
            snpeffann_entries = snpeffann.split(',')

            # For the moment, it only takes the first entry
            # Get the effect
            effect, first_snpeff_entry = snpeffann_entries[0].split('(')

            # first_snpeff_entry, _ = first_snpeff_entry.split(')')
            #
            # Separete each item from the first entry
            # It might change if I implement another role for effects
            # This is not being used for the simple rule of taking
            # the first entry.
            # Effect_Impact, Functional_Class, Codon_Change, Amino_Acid_Change,
            # Amino_Acid_length, Gene_Name, Transcript_BioType, Gene_Coding,
            # Transcript_ID, Exon_Rank, Genotype, *errors = first_snpeff_entry.split('|')
            
            # Re-assemble the EFF field with the simplified version
            # of the SNPEff annottion.
            fields[7] = altcount + ";" + refcount + ";" + "EFF=" + snpeffann_entries[0]

            # Re-assemble the entire line
            reassembled_line = '\t'.join(fields)

            if keeponlyterms:
                if any(x == effect for x in relevant_effect_terms):
                    output_file.write(reassembled_line + "\n")
            else:
                output_file.write(reassembled_line + "\n")

    return "file processed"


def parseargs():
    parser = argparse.ArgumentParser("python simplify_snpeff.py", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", help="The the name of the file which the data are to be read from (annotated vcf from SNPEff)",
                        dest="file", required=True, type=str)
    parser.add_argument("-o", help="A character string naming a file",
                        dest="outfile", default=None, type=str)
    parser.add_argument("-f", help="Keep only relevant terms (it removes up(down)-stream SNPs)",
                        dest="keepterms", default=False, type=bool)
    parser.add_argument("-m", help="Method to consolidate SNPEff annotation",
                        dest="method", default="First", type=str)
    return parser


def main(argv):
    """
    This is the main program definition.
    """
    parser = parseargs()
    if argv[-1] == '':
        argv = argv[0:-1]
    args = parser.parse_args(argv)

    file = args.file
    outfile = args.outfile
    keepterms = args.keepterms
    method = args.method

    result = simplify_snpeff(
        file=file, outfile=outfile, keeponlyterms=keepterms, method=method
    )
    print(result)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        main(['-h'])
    else:
        main(sys.argv[1:])












# Majority-rule like effects filtering
# Create an empty dictionary to store the counts
# item_counts = {}

# # Loop through the list and count the occurrences of each item
# for item in snpeffann_effects_entries:
#     if item in item_counts:
#         item_counts[item] += 1
#     else:
#         item_counts[item] = 1

# # Use a threshold to define intron_variant to save
# # Calculate the total count of all items
# total_count = sum(item_counts.values())

# # Calculate the total count of the items you want to retain
# total_retain_count = sum(item_counts.get(item, 0) for item in items_to_retain)

# # Calculate the percentage for each item you want to retain
# item_percentages = {item: (item_counts.get(item, 0) / total_count) * 100 for item in items_to_retain}

# # Check if each item's percentage is >= 30%
# valid_items = [item for item in items_to_retain if item_percentages.get(item, 0) >= 30]

# # If all items meet the condition, retain them; otherwise, don't retain any
# if len(valid_items) > 0:
# #print(line)
# filtered_vcf_lines.append(line)

