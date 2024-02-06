"""
Program converts an annotated and simplified vcf file to a .tsv table
"""
import os
import subprocess
import argparse
import sys
from typing import Tuple


# Count the number of alternative and reference genotypes
def count_genotypes(genotypes: list) -> Tuple[int, int]:
    """
    Count the number of alternative and reference genotypes
    in a list of genotypes and return a string with the counts
    """
    refcount = 0
    altcount = 0
    for genotype in genotypes:
        if genotype != './.':
            if genotype == "0/0":
                refcount += 1
            elif genotype == "1/1":
                altcount += 1

    return altcount, refcount


# Process the first line of the input file
def process_first_line(inputfile: str) -> Tuple[str, str, str]:
    """
    Process the first line of the input file. This is used in conjunction
    with a railroad pattern to check if the input file has the right format.
    """
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

        # Return the results for the railroad pattern
        return first, second, third
    else:
        raise ValueError("Input file is empty")


# Process refaltcount VCF
def processes_refaltcount_vcf(inputfile: str, chrm_list: list,
                              reference: str, samtools: str, nflankinbps: int
                              ) -> Tuple[list, list, list]:
    """
    Railroad pattern #1: INFO field has ALTCOUNT AND REFCOUNT.
    """

    print("Start processing the annotate VCF file. It might take a while...")
    with open(inputfile, "r", encoding="utf-8") as input_file:
        current_block = []
        lines_in_block = []
        list_block = []
        list_lines_in_block = []
        previous_position = 0
        line_trck = 0
        # Read each line in the file and save to a list
        lines = []
        for line in input_file:
            if line.startswith("##") or line.startswith("#"):
                continue

            # Extract variant information from the VCF fields
            fields = line.strip().split("\t")

            if len(fields) < 8:
                continue

            # Unpack FIELDS items
            # anc/ref, derived/alt
            chrom, pos, vcfid, ref, alt, qual, filtr, info = fields[
                :8
            ]  # maybe use pop() here

            # Cast pos to integer
            pos = int(pos)

            # Unpack items in info
            altcount_, refcount_, snpeff_, *sift_ = info.strip().split(";")

            _, altcount = altcount_.strip().split("=")
            _, refcount = refcount_.strip().split("=")

            # Get the mutational context
            # Run samtools faidx as a python subprocess
            refflank, altflank = "NA", "NA"
            refcodon, altcodon, refaa, altaa = "NA", "NA", "NA", "NA"
            effect = "NA"
            (
                effect_impact,
                functional_class,
                gene_name,
                transcript_biotype,
                gene_coding,
                transcript_id,
            ) = ("NA", "NA", "NA", "NA", "NA", "NA")
            (
                transcript,
                geneid,
                genename,
                region,
                varianttype,
                siftscore,
                siftmedian,
                siftpred,
                deleteriousness,
            ) = ("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")

            # Frim this part until the end, run only in canonical chromosomes
            if any(x == chrom for x in chrm_list):
                nt_bf = subprocess.run(
                    [
                        samtools,
                        "faidx",
                        reference,
                        (chrom + ":" + str(pos - nflankinbps) + "-" + str(pos - 1)),
                    ],
                    capture_output=True, check=False
                )
                nt_af = subprocess.run(
                    [
                        samtools,
                        "faidx",
                        reference,
                        (chrom + ":" + str(pos + 1) + "-" + str(pos + nflankinbps)),
                    ],
                    capture_output=True, check=False
                )

                _, flkng_bf, *_ = str(nt_bf.stdout).strip().split("\\n")
                _, flkng_af, *_ = str(nt_af.stdout).strip().split("\\n")

                refflank = flkng_bf.upper() + ref + flkng_af.upper()
                altflank = flkng_bf.upper() + alt + flkng_af.upper()

                # This part is for fixing basis in flaking bases string
                # if SNPs are close (less than nflankinbps)
                # This keeps track of block of positions with near SNPs
                current_position = pos

                if (
                    current_position - previous_position
                ) < nflankinbps:  # (maybe nflankinbps - 1)
                    if not current_block:
                        current_block.append(previous_position)
                        lines_in_block.append(line_trck - 1)
                    current_block.append(current_position)
                    lines_in_block.append(line_trck)
                else:
                    current_block = []
                    lines_in_block = []

                if len(current_block) > 1:
                    if current_block not in list_block:
                        list_block.append(current_block)
                    if lines_in_block not in list_lines_in_block:
                        list_lines_in_block.append(lines_in_block)

                previous_position = current_position

                # Unpack items in SNPeff info::eff
                _, eff_ = snpeff_.strip().split("=")
                effect, eff = eff_.strip().split("(")  # (gatk format)
                eff, _ = eff.strip().split(")")

                # ALL SNPEff fields (gatk format)
                (
                    effect_impact,
                    functional_class,
                    codon_change,
                    aa_change,
                    aa_length,
                    gene_name,
                    transcript_biotype,
                    gene_coding,
                    transcript_id,
                    exon,
                    genotype,
                    *errors,
                ) = eff.strip().split("|")

                # Unpack codons (for gatk format)
                # if effect == "SYNONYMOUS_CODING" or effect == "NON_SYNONYMOUS_CODING":
                if effect in ("SYNONYMOUS_CODING", "NON_SYNONYMOUS_CODING"):
                    refcodon, altcodon = codon_change.strip().split("/")

                # Unpack items in info::sift (when present)
                if len(sift_) > 0:
                    # Unpack items in info::eff
                    _, siftinfo = sift_[0].strip().split("=")

                    # ALL SIFT fields
                    (
                        allele,
                        transcript,
                        geneid,
                        genename,
                        region,
                        varianttype,
                        aa,
                        aaposition,
                        siftscore,
                        siftmedian,
                        siftnumseqs,
                        alleletype,
                        siftpred,
                    ) = siftinfo.strip().split("|")
                    refaa, altaa = aa.strip().split("/")

                    # Define deleteriousness status
                    if effect == "NON_SYNONYMOUS_CODING":
                        if siftscore == "NA":
                            deleteriousness = "NA"
                        else:
                            if float(siftscore) < 0.05:
                                deleteriousness = "deleterious"
                            else:
                                deleteriousness = "tolerated"
                    elif effect == "SYNONYMOUS_CODING":
                        deleteriousness = "NA"

            # Create a re-usable list to store the needed information
            line_list = [
                chrom,
                pos,
                vcfid,
                ref,
                alt,
                refcount,
                altcount,
                refflank,
                altflank,
                refcodon,
                altcodon,
                refaa,
                altaa,
                effect,
                effect_impact,
                functional_class,
                gene_name,
                transcript_biotype,
                gene_coding,
                transcript_id,
                transcript,
                geneid,
                genename,
                region,
                varianttype,
                siftscore,
                siftmedian,
                siftpred,
                deleteriousness,
            ]
            lines.append(line_list)
            line_trck += 1

    return lines, list_lines_in_block, list_block


# Process the standard VCF
def processes_standard_vcf(inputfile: str, chrm_list: list,
                           reference: str, samtools: str, nflankinbps: int
                           ) -> Tuple[list, list, list]:
    """
    Railroad pattern #2: INFO field has AC AND AF.
    Genotype counts need to be calculated with count_genotypes().
    """

    print("Start processing the annotate VCF file. It might take a while...")
    with open(inputfile, "r", encoding="utf-8") as input_file:
        current_block = []
        lines_in_block = []
        list_block = []
        list_lines_in_block = []
        previous_position = 0
        line_trck = 0
        # Read each line in the file and save to a list
        lines = []
        for line in input_file:
            if line.startswith("##") or line.startswith("#"):
                continue

            # Extract variant information from the VCF fields
            fields = line.strip().split("\t")

            if len(fields) < 8:
                continue

            # Unpack FIELDS items
            # anc/ref, derived/alt
            chrom, pos, vcfid, ref, alt, qual, filtr, info = fields[
                :8
            ]  # maybe use pop() here

            genotypes = fields[9:]

            # Cast pos to integer
            pos = int(pos)

            # Unpack items in info
            ac_, af_, snpeff_, *sift_ = info.strip().split(";")

            # Get the mutational context
            # Run samtools faidx as a python subprocess
            refflank, altflank = "NA", "NA"
            refcodon, altcodon, refaa, altaa = "NA", "NA", "NA", "NA"
            effect = "NA"
            (
                effect_impact,
                functional_class,
                gene_name,
                transcript_biotype,
                gene_coding,
                transcript_id,
            ) = ("NA", "NA", "NA", "NA", "NA", "NA")
            (
                transcript,
                geneid,
                genename,
                region,
                varianttype,
                siftscore,
                siftmedian,
                siftpred,
                deleteriousness,
            ) = ("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")

            # Frim this part until the end, run only in canonical chromosomes
            if any(x == chrom for x in chrm_list):
                nt_bf = subprocess.run(
                    [
                        samtools,
                        "faidx",
                        reference,
                        (chrom + ":" + str(pos - nflankinbps) + "-" + str(pos - 1)),
                    ],
                    capture_output=True, check=False
                )
                nt_af = subprocess.run(
                    [
                        samtools,
                        "faidx",
                        reference,
                        (chrom + ":" + str(pos + 1) + "-" + str(pos + nflankinbps)),
                    ],
                    capture_output=True, check=False
                )

                _, flkng_bf, *_ = str(nt_bf.stdout).strip().split("\\n")
                _, flkng_af, *_ = str(nt_af.stdout).strip().split("\\n")

                refflank = flkng_bf.upper() + ref + flkng_af.upper()
                altflank = flkng_bf.upper() + alt + flkng_af.upper()

                # This part is for fixing basis in flaking bases string
                # if SNPs are close (less than nflankinbps)
                # This keeps track of block of positions with near SNPs
                current_position = pos

                if (
                    current_position - previous_position
                ) < nflankinbps:  # (maybe nflankinbps - 1)
                    if not current_block:
                        current_block.append(previous_position)
                        lines_in_block.append(line_trck - 1)
                    current_block.append(current_position)
                    lines_in_block.append(line_trck)
                else:
                    current_block = []
                    lines_in_block = []

                if len(current_block) > 1:
                    if current_block not in list_block:
                        list_block.append(current_block)
                    if lines_in_block not in list_lines_in_block:
                        list_lines_in_block.append(lines_in_block)

                previous_position = current_position

                # Unpack items in SNPeff info::eff
                _, eff_ = snpeff_.strip().split("=")
                effect, eff = eff_.strip().split("(")  # (gatk format)
                eff, _ = eff.strip().split(")")

                # ALL SNPEff fields (gatk format)
                (
                    effect_impact,
                    functional_class,
                    codon_change,
                    aa_change,
                    aa_length,
                    gene_name,
                    transcript_biotype,
                    gene_coding,
                    transcript_id,
                    exon,
                    genotype,
                    *errors,
                ) = eff.strip().split("|")

                # Unpack codons (for gatk format)
                # if effect == "SYNONYMOUS_CODING" or effect == "NON_SYNONYMOUS_CODING":
                if effect in ("SYNONYMOUS_CODING", "NON_SYNONYMOUS_CODING"):
                    refcodon, altcodon = codon_change.strip().split("/")

                # Unpack items in info::sift (when present)
                if len(sift_) > 0:
                    # Unpack items in info::eff
                    _, siftinfo = sift_[0].strip().split("=")

                    # ALL SIFT fields
                    (
                        allele,
                        transcript,
                        geneid,
                        genename,
                        region,
                        varianttype,
                        aa,
                        aaposition,
                        siftscore,
                        siftmedian,
                        siftnumseqs,
                        alleletype,
                        siftpred,
                    ) = siftinfo.strip().split("|")
                    refaa, altaa = aa.strip().split("/")

                    # Define deleteriousness status
                    if effect == "NON_SYNONYMOUS_CODING":
                        if siftscore == "NA":
                            deleteriousness = "NA"
                        else:
                            if float(siftscore) < 0.05:
                                deleteriousness = "deleterious"
                            else:
                                deleteriousness = "tolerated"
                    elif effect == "SYNONYMOUS_CODING":
                        deleteriousness = "NA"

            altcount, refcount = count_genotypes(genotypes)
            
            # Create a re-usable list to store the needed information
            line_list = [
                chrom,
                pos,
                vcfid,
                ref,
                alt,
                refcount,
                altcount,
                refflank,
                altflank,
                refcodon,
                altcodon,
                refaa,
                altaa,
                effect,
                effect_impact,
                functional_class,
                gene_name,
                transcript_biotype,
                gene_coding,
                transcript_id,
                transcript,
                geneid,
                genename,
                region,
                varianttype,
                siftscore,
                siftmedian,
                siftpred,
                deleteriousness,
            ]
            lines.append(line_list)
            line_trck += 1

    return lines, list_lines_in_block, list_block

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

    # Check if input files 
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
                current_refflanking = list(lines[current_element][7])
                current_altflanking = list(lines[current_element][8])

                # Check if there are multiple elements before or after
                if len(elements_before) >= 1:
                    for j, eb in enumerate(elements_before):
                        dist = abs(current_pos - pos_before[j])
                        if dist > nflankinbps:
                            continue
                        else:
                            ref_before = lines[eb][3]
                            alt_before = lines[eb][4]
                            current_refflanking[nflankinbps - dist] = ref_before
                            current_altflanking[nflankinbps - dist] = alt_before

                if len(elements_after) >= 1:
                    for j, ea in enumerate(elements_after):
                        dist = abs(current_pos - pos_after[j])
                        print(dist)
                        if dist > nflankinbps:
                            continue
                        else:
                            ref_after = lines[ea][3]
                            alt_after = lines[ea][4]
                            current_refflanking[nflankinbps + dist] = ref_after
                            current_altflanking[nflankinbps + dist] = alt_after

                # Here update the flanking sequences in the corresponding line
                lines[current_element][7] = "".join(current_refflanking)
                lines[current_element][8] = "".join(current_altflanking)

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
