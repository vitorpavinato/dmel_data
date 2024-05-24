"""
Program converts an annotated and simplified vcf file to a .tsv table
"""

import argparse
import sys
from vcf_to_tsv_utils import *
from vcf_to_tsv_implementations import *


inputfile = "/Users/tur92196/WorkDir/DGN/dgrp2/masked/vcfs/sift4g/NC_Chr2L_sample.vcf"
outputfile = "annotate_vcf/examples/tables/example_remade_rooted_lifted_ann_simplified_SIFTpredictions_table.vcf"
samtools_path = "/Users/tur92196/local/samtools1_8/bin/samtools"
reference = "/Users/tur92196/WorkDir/DGN/reference/dm6.fa"

nflankinbps = 4
custom_effect_name = "dm6_short_introns"
new_custom_effect_name = "SI"
sift4g_annotations = True

# vcf_to_tsv.py railroad pattern implementation
def vcf_to_tsv(
    inputfile: str, outputfile: str, 
    reference: str, samtools_path: str, nflankinbps: int,
    custom_effect_name: str = None, new_custom_effect_name: str = None, 
    sift4g_annotations: bool = False
) -> None:
    """Scratching the code..."""

    # Input file check
    inputfile = check_input_file(inputfile)

    # Output file check
    outputfile = output_file_name(inputfile, outputfile)

    # Samtools path check
    samtools = check_samtools_path(samtools_path)

    # Reference file check
    reference = check_reference_genome_file(reference)

    # Create .fai file if it doesn't exist
    create_faidx(reference, samtools)

    # Here is the railroad pattern implementation
    if sift4g_annotations:
        vcf_lines, list_lines_in_block, list_block = processes_snpeff_sift4g_vcf(
            inputfile=inputfile, reference=reference,
            samtools=samtools, nflankinbps=nflankinbps,
            custom_effect_name=custom_effect_name,
            new_custom_effect_name=new_custom_effect_name
        )
        header = snpeff_sift4g_header()

    else:
        vcf_lines, list_lines_in_block, list_block = processes_snpeff_vcf(
            inputfile=inputfile, reference=reference,
            samtools=samtools, nflankinbps=nflankinbps,
            custom_effect_name=custom_effect_name,
            new_custom_effect_name=new_custom_effect_name
        )
        header = snpeff_header()

    # Execute the rest of the function.
    print("Fixing ref and alt mutations on flanking bases. It might take a while...")
    # for block_lines, block_pos in zip(new_list_lines_in_block, new_list_block):
    for block_lines, block_pos in zip(list_lines_in_block, list_block):
        blocksize = len(block_lines)
        if blocksize > 1:
            for i, current_element in enumerate(block_lines):
                elements_before = block_lines[:i]
                elements_after = block_lines[i + 1:]
                current_pos = block_pos[i]
                pos_before = block_pos[:i]
                pos_after = block_pos[i + 1:]

                # Get the flanking bases that need to be updated
                current_refcontext = list(vcf_lines[current_element][13]) # refcontext
                current_altcontext = list(vcf_lines[current_element][14]) # altcontext

                # Check if there are multiple elements before or after
                if len(elements_before) >= 1:
                    for j, eb in enumerate(elements_before):
                        dist = abs(current_pos - pos_before[j])
                        if dist > nflankinbps:
                            continue
                        else:
                            ref_before = vcf_lines[eb][3] # ref before
                            alt_before = vcf_lines[eb][4] # alt before
                            current_refcontext[nflankinbps - dist] = ref_before
                            current_altcontext[nflankinbps - dist] = alt_before

                if len(elements_after) >= 1:
                    for j, ea in enumerate(elements_after):
                        dist = abs(current_pos - pos_after[j])
                        print(dist)
                        if dist > nflankinbps:
                            continue
                        else:
                            ref_after = vcf_lines[ea][3] # ref after
                            alt_after = vcf_lines[ea][4] # alt after
                            current_refcontext[nflankinbps + dist] = ref_after
                            current_altcontext[nflankinbps + dist] = alt_after

                # Here update the flanking sequences in the corresponding line
                vcf_lines[current_element][13] = "".join(current_refcontext)
                vcf_lines[current_element][14] = "".join(current_altcontext)
    
    
    # Write .tsv file
    write_tsv_file(vcf_lines, header, outputfile)





# This function converts entries on a VCF file to rows of a table
# This is a higher level function
def vcf_to_tsv(file: str = None, outfile: str = None,
               reference: str = None,
               nflankinbps: int = 4,
               samtools_path: str = None) -> None:
    """
    This function converts entries on a VCF file to rows of a table.
    It annotates the mutational context of SNPs of Drosophia melanogaster.
    It expects the same reference file used in the VCF (downloaded or lifted over).
    """

    # Check if input files are provided
    if file == "" or file is None:
        raise ValueError("file must be provided")

    # Check if the file exists
    if not os.path.exists(file):
        raise ValueError("file does not exist")

    # Check if samtools path is provided
    if samtools_path is None:
        raise ValueError("samtools_path must be provided")

    samtools = samtools_path

    # Check if the reference file is provided
    if reference is None:
        raise ValueError("reference must be provided")

    # Check if faidx associated with the reference exists
    if not os.path.exists((reference + ".fai")):
        faidx = subprocess.run([samtools, "faidx", reference], capture_output=True, check=False)

    # Get the path of the input file
    path, filename = os.path.split(file)

    # Define the basename for the the output file in case it is not provided
    basename, _ = filename.strip().split("_simplified_SIFTpredictions.vcf")

    # Define input and output files
    inputfile = file

    if outfile is None:
        outputfile = path + "/" + basename + "_table.tsv"
    else:
        outputfile = outfile

    # Define the list of chromosomes and the path to the reference file based 
    # on the selected species
    chrm_list = ["chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"]

    # This defines the header of the output file
    header = [
        "chrom",
        "pos",
        "id",
        "ref",
        "alt",
        "refcount",
        "altcount",
        "refflank",
        "altflank",
        "refcodon",
        "altcodon",
        "refaa",
        "altaa",
        "effect",
        "snpeff_effimpact",
        "snpeff_funclass",
        "snpeff_genename",
        "snpeff_trnscbiotype",
        "snpeff_genecoding",
        "snpeff_trnscid",
        "sift_trnscid",
        "sift_geneid",
        "sift_genename",
        "sift_region",
        "sift_vartype",
        "sifts_core",
        "sift_median",
        "sift_pred",
        "deleteriousness",
    ]
    header = "\t".join(str(item) for item in header)

    # Implemeting railroads pattern design
    first, second, third = process_first_line(inputfile)

    if first == "ALTCOUNT" and second == "REFCOUNT":
        # Execute code for the ALTCOUNT and REFCOUNT condition
        lines, list_lines_in_block, list_block = processes_refaltcount_vcf(inputfile, chrm_list, reference, samtools, nflankinbps)

    elif first == "AC" and second == "AF":
        # Execute code for the AC and AF condition
        lines, list_lines_in_block, list_block = processes_standard_vcf(inputfile, chrm_list, reference, samtools, nflankinbps)
    else:
        # Handle the case where the first and second values do not match the expected conditions
        print("Error: Unsupported values for first and second fields")
        raise ValueError("Unsupported values for first and second fields")

    # Execute the rest of the function.
    print("Fixing ref and alt mutations on flanking bases. It might take a while...")
    # for block_lines, block_pos in zip(new_list_lines_in_block, new_list_block):
    for block_lines, block_pos in zip(list_lines_in_block, list_block):
        blocksize = len(block_lines)
        if blocksize > 1:
            for i, current_element in enumerate(block_lines):
                elements_before = block_lines[:i]
                elements_after = block_lines[i + 1:]
                current_pos = block_pos[i]
                pos_before = block_pos[:i]
                pos_after = block_pos[i + 1:]

                # Get the flanking bases that need to be updated
                current_refcontext = list(lines[current_element][7])
                current_altcontext = list(lines[current_element][8])

                # Check if there are multiple elements before or after
                if len(elements_before) >= 1:
                    for j, eb in enumerate(elements_before):
                        dist = abs(current_pos - pos_before[j])
                        if dist > nflankinbps:
                            continue
                        else:
                            ref_before = lines[eb][3]
                            alt_before = lines[eb][4]
                            current_refcontext[nflankinbps - dist] = ref_before
                            current_altcontext[nflankinbps - dist] = alt_before

                if len(elements_after) >= 1:
                    for j, ea in enumerate(elements_after):
                        dist = abs(current_pos - pos_after[j])
                        print(dist)
                        if dist > nflankinbps:
                            continue
                        else:
                            ref_after = lines[ea][3]
                            alt_after = lines[ea][4]
                            current_refcontext[nflankinbps + dist] = ref_after
                            current_altcontext[nflankinbps + dist] = alt_after

                # Here update the flanking sequences in the corresponding line
                lines[current_element][7] = "".join(current_refcontext)
                lines[current_element][8] = "".join(current_altcontext)

    # Prompt message
    print("Exporting the processed VCF to a .TSV file")
    with open(outputfile, "a", encoding="utf-8") as fo:
        fo.write(header + "\n")
        for line in lines:
            fo.write(("\t".join(str(item) for item in line)) + "\n")
    # print("DONE!!!")

    return "file processed"


# Define the command line arguments
def parseargs():
    """
    Function defines command-line parsing arguments.
    """
    parser = argparse.ArgumentParser(
        "python vcf_to_table.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-f",
        help="The string name of the input vcf file (with annotation; with the path to)",
        dest="file",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-o",
        help="The string name of the tsv output file (with the path to)",
        dest="outfile",
        default=None,
        type=None,
    )
    parser.add_argument(
        "-r",
        help="Path to the reference genome of the vcf file",
        dest="reference",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-n",
        help="Number of bases flanking each targeted SNP",
        dest="nflankinbps",
        default=4,
        type=int,
    )
    parser.add_argument(
        "-s",
        help="Path to the samtools",
        dest="samtools_path",
        required=True,
        type=str,
    )
    return parser


def main(argv) -> None:
    """
    This is the main program definition.
    """
    parser = parseargs()
    if argv[-1] == "":
        argv = argv[0:-1]
    args = parser.parse_args(argv)

    # Define input and output files
    file = args.file
    outfile = args.outfile
    reference = args.reference
    nflankinbps = args.nflankinbps
    samtools_path = args.samtools_path

    result = vcf_to_tsv(file=file, outfile=outfile, reference=reference, nflankinbps=nflankinbps, samtools_path=samtools_path)
    print(result)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        main(["-h"])
    else:
        main(sys.argv[1:])
