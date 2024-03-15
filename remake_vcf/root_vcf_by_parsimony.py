"""
Annotate SNPs with the ancestral allele by parsimony.
"""

import sys
import argparse
import subprocess  # this is being used to run samtools
from pathlib import Path  # This is being used to check if the file exists
from typing import List
from remake_vcf import create_samtools_faix


def root_snps_by_parsimony(
    input_file: Path = None,
    output_file: Path = None,
    outgroup_file: Path = None,
    samtools_path: Path = None
) -> None:
    """
    Root SNPs in a vcf based on the outgroup sequence.
    It takes a vcf file and return a new vcf with the rooted SNPs.
    """

    # Check if all required arguments are provided
    # This will raise an error if any of the required arguments are None

    # This is for the input_file
    if input_file is None:
        raise ValueError("input_file must be provided")

    # This is for the output_file
    if output_file is None:
        raise ValueError("output_file must be provided")

    if outgroup_file is None:
        raise ValueError("outgroup_file must be provided")

    # This is for the samtools_path
    if samtools_path is None:
        raise ValueError("samtools_path must be provided")

    # Check the existence of a outgroup reference file faidx
    create_samtools_faix(outgroup_file, samtools_path)
    print("Faidx created for the outgroup reference...")

    # Copy the input vcf header to the output vcf
    with open(input_file, "r", encoding="utf-8") as input_file_obj, open(output_file, "a", encoding="utf-8") as output_file_obj:
        for line in input_file_obj:
            if (line.startswith("##") or line.startswith("#")):
                output_file_obj.write(line)
            else:
                break

    # Process remaining lines
    with open(input_file, "r", encoding="utf-8") as input_file_obj, open(output_file, "a", encoding="utf-8") as output_file_obj:
        for line in input_file_obj:
            if (line.startswith("##") or line.startswith("#")):
                continue

            # Extract variant information from the VCF fields
            fields = line.strip().split("\t")

            # Skip empty lines or lines with less than 9 fields
            if len(fields) < 9:
                continue

            # Create a list of alleles
            alleles = []

            # First add the reference allele
            alleles.append(fields[3])

            # Then extend the list to include the alternative allele(s)
            alleles.extend(fields[4].split(","))

            # Call samtools faidx to check outgroup reference base
            outgroup_ref_base = subprocess.run([
                samtools_path,
                "faidx",
                outgroup_file,
                (fields[0] + ":" + fields[1] + "-" + fields[1]),
            ],
                                    capture_output=True,
                                    check=False)

            # Get the reference base as a variable
            _, outgroup_ref_base, *_ = str(outgroup_ref_base.stdout).strip().split("\\n")

            # This will make sure that masked genome will work
            outgroup_ref_base = outgroup_ref_base.upper()

            if outgroup_ref_base in alleles:
                fields[2] = f"AA={outgroup_ref_base}"
            else:
                fields[2] = f"AA={fields[3]}"

            # Write the new line
            reassembled_line = '\t'.join(fields)

            # Write the new line
            with open(output_file, 'a', encoding='utf-8') as output_file_obj:
                output_file_obj.write(reassembled_line + "\n")


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser("python root_vcf_by_parsimony.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_file", type=str, required=True, help="Input vcf file", dest="input_file")
    parser.add_argument("-o", "--output_file", type=str, required=True, help="Output vcf file", dest="output_file")
    parser.add_argument("-r", "--outgroup_ref_file", type=str, required=True, help="Outgroup reference file", dest="outgroup_file")
    parser.add_argument("-s", "--samtools_path", type=str, required=True, help="Path to samtools", dest="samtools_path")
    return parser


def main(argv: List[str]) -> None:
    """
    Main function
    """
    # Parse command line arguments
    parser = parse_args()
    if argv[-1] == '':
        argv = argv[0:-1]
    args = parser.parse_args(argv)

    # Define input and output files
    input_file = Path(args.input_file)
    output_file = Path(args.output_file)
    outgroup_file = Path(args.outgroup_file)
    samtools_path = Path(args.samtools_path)

    # Run root_vcf_by_parsimony
    root_snps_by_parsimony(input_file, output_file, outgroup_file, samtools_path)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        main(["-h"])
    else:
        main(sys.argv[1:])
