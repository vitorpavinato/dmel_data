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
def simplify_snpeff(file: str, outfile: str = None,
                    keeponlyterms: bool = False,
                    method: str = "First") -> None:

    """
    This function simplifies SNPEff entries of a VCF file.
    For the moment, method don't do anything.
    """

    path, filename = os.path.split(file)

    # Define the basename for the outputs from the filename
    basename, _ = filename.strip().split('.vcf')

    # Define input and output files
    inputfile = file

    # Check if the input file has the right INFO fields format
    with open(inputfile, "r", encoding="utf-8") as input_file:
        first_line = None
        for line in input_file:
            if not line.startswith("##") and not line.startswith("#"):
                first_line = line
                break

    if first_line is not None:
        fields = first_line.strip().split('\t')
        # Check the INFO fields
        info_field = fields[7]
        info_field_list = info_field.split(';')

        # Check if the INFO fields have at least 3 elements or if it is not empty
        if len(info_field_list) < 3 or info_field == "":
            raise ValueError("Input file has missing elements in INFO field. Was this vcf annotated with SNPEff?")

        # Handle the case when the file is not empty
        # and has at least 3 elements in INFO field
        # Split the elements in info_field_list
        first = info_field_list[0].split('=')[0]
        second = info_field_list[1].split('=')[0]
        third = info_field_list[2].split('=')[0]

        if (first != "AC" or first != "ALTCOUNT") and (second != "AF" or second != "REFCOUNT") and (third != "EFF" ):
            raise ValueError("Input is not supported by this script!")

    else:
        # Handle the case when the file is empty
        raise ValueError("Input file has only header information...")

    # This handles the output file name
    if outfile is None:
        outputfile = path + '/' + basename + "_simplified.vcf"
    else:
        outputfile = outfile

    # This handles if keeponlyterms is set to True
    if keeponlyterms:
        # Define annotation effects terms we care about
        relevant_effect_terms = ['INTRON', 'SYNONYMOUS_CODING', 'NON_SYNONYMOUS_CODING', 'INTERGENIC']

    # Copy the header from the inputfile to the ouputfile
    with open(inputfile, "r", encoding="utf-8") as input_file, open(outputfile, "w", encoding="utf-8") as output_file:
        for line in input_file:
            if (line.startswith("##") or line.startswith("#")):
                output_file.write(line)

    # Open the input file in read mode and output file in write mode
    with open(inputfile, "r") as input_file, open(outputfile, "a") as output_file:
        # For each line, breakdow the SNPEff info field to retain only the effect entries
        # Iterate through the lines in the input file
        for line in input_file:
            if (line.startswith("##") or line.startswith("#")):
                continue
        
            # Split the line into fields using tab as the delimiter
            fields = line.strip().split('\t')

            # Get the ANN field (assuming it's always the eighth INFO field)
            info_field = fields[7]

            # Split the ANN field by commas to separete the SNPEff annotation
            eli, elj, snpeffs, *lo = info_field.split(';')

            # Check if any ReverseComplementedAlleles is present
            # in any SNPEff annotation entry and deal with it
            if (snpeffs == "ReverseComplementedAlleles"):
                _, snpeffann, = lo[0].split("EFF=")
            else:
                _, snpeffann, = snpeffs.split("EFF=")

            # Separete each SNPEff entry
            snpeffann_entries = snpeffann.split(',')

            # It only takes the first entry
            # Get the effect to use latter in keeping only relevant terms
            effect, first_snpeff_entry = snpeffann_entries[0].split('(')
            
            # Re-assemble the EFF field with the simplified version
            # of the SNPEff annottion.
            fields[7] = eli + ";" + elj + ";" + "EFF=" + snpeffann_entries[0]

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
                        dest="keeponlyterms", default=False, type=bool)
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
    keeponlyterms = args.keeponlyterms
    method = args.method

    result = simplify_snpeff(
        file=file, outfile=outfile, keeponlyterms=keeponlyterms, method=method
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

