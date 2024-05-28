"""
Program converts an annotated and simplified vcf file to a .tsv table
"""

import argparse
import sys
from vcf_to_tsv_func import *
from implementations import *
from mutational_context_func import fix_mutational_context


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
    vcf_lines = fix_mutational_context(
        list_lines_in_block=list_lines_in_block, 
        list_block=list_block,
        vcf_lines=vcf_lines, 
        nflankinbps=nflankinbps
    )

    # Write .tsv file
    write_tsv_file(vcf_lines, header, outputfile)

    return "file processed"


def main():
    """
    SOMETHING HERE
    """
    
    # This was to test in larger file which might be easy to find closest SNPs
    # inputfile = "/Users/tur92196/WorkDir/DGN/dgrp2/masked/vcfs/sift4g/NC_Chr2L_sample.vcf"
    inputfile = "annotate_vcf/examples/sift4g/example_remade_rooted_lifted_ann_simplified_SIFTpredictions.vcf"
    outputfile = "annotate_vcf/examples/tables/example_remade_rooted_lifted_ann_simplified_SIFTpredictions_table.vcf"
    samtools_path = "/Users/tur92196/local/samtools1_8/bin/samtools"
    reference = "/Users/tur92196/WorkDir/DGN/reference/dm6.fa"

    nflankinbps = 4
    custom_effect_name = "dm6_short_introns"
    new_custom_effect_name = "SI"
    sift4g_annotations = True

    result = vcf_to_tsv(inputfile=inputfile, outputfile=outputfile,
                        reference=reference, samtools_path=samtools_path,
                        nflankinbps=nflankinbps, custom_effect_name=custom_effect_name,
                        new_custom_effect_name=new_custom_effect_name,
                        sift4g_annotations=sift4g_annotations)

    print(result)


if __name__ == "__main__":
    main()




# # Define the command line arguments
# def parseargs():
#     """
#     Function defines command-line parsing arguments.
#     """
#     parser = argparse.ArgumentParser(
#         "python vcf_to_table.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter
#     )
#     parser.add_argument(
#         "-f",
#         help="The string name of the input vcf file (with annotation; with the path to)",
#         dest="file",
#         required=True,
#         type=str,
#     )
#     parser.add_argument(
#         "-o",
#         help="The string name of the tsv output file (with the path to)",
#         dest="outfile",
#         default=None,
#         type=None,
#     )
#     parser.add_argument(
#         "-r",
#         help="Path to the reference genome of the vcf file",
#         dest="reference",
#         required=True,
#         type=str,
#     )
#     parser.add_argument(
#         "-n",
#         help="Number of bases flanking each targeted SNP",
#         dest="nflankinbps",
#         default=4,
#         type=int,
#     )
#     parser.add_argument(
#         "-s",
#         help="Path to the samtools",
#         dest="samtools_path",
#         required=True,
#         type=str,
#     )
#     return parser


# def main(argv) -> None:
#     """
#     This is the main program definition.
#     """
#     parser = parseargs()
#     if argv[-1] == "":
#         argv = argv[0:-1]
#     args = parser.parse_args(argv)

#     # Define input and output files
#     file = args.file
#     outfile = args.outfile
#     reference = args.reference
#     nflankinbps = args.nflankinbps
#     samtools_path = args.samtools_path

#     result = vcf_to_tsv(file=file, outfile=outfile, reference=reference, nflankinbps=nflankinbps, samtools_path=samtools_path)
#     print(result)


# if __name__ == "__main__":
#     if len(sys.argv) < 2:
#         main(["-h"])
#     else:
#         main(sys.argv[1:])
