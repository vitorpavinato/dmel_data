"""
This program converts an annotated and simplified
vcf file to a .tsv table
"""
import os
import subprocess
import argparse
import sys
from enum import Enum


# CLASSES
class Species(Enum):
    """
    Class to create a enum for a species list
    """
    DMEL = 1
    HSAP = 2


# FUNCTIONS
def vcf_to_tsv(file: str, outfile: str,
               species: str = "Dmel", nflankinbps: int = 4) -> None:
    """
    This function converts entries on a VCF file to rows of a table.
    Because it annotates the mutational context of SNPs, for the moment it
    only works for two model organims, Drosophia melanogaster and Humans.
    Maybe add an option for a user specified reference file and annotation.
    """

    # Config
    samtools = "/Users/tur92196/local/samtools1_8/bin/samtools"

    path, filename = os.path.split(file)

    # Define the basename for the outputs from the filename
    basename, _ = filename.strip().split("_simplified_SIFTpredictions_")

    # Define input and output files
    inputfile = file

    if outfile is None:
        outputfile = path + basename + "_table.tsv"
    else:
        outputfile = outfile

    # Define the list of chromosomes and the path to the reference file based 
    # on the selected species
    if species == Species.DMEL.name:
        chrm_list = ["2L", "2R", "3L", "3R", "4", "X", "Y"]
        reference = "/Users/tur92196/WorkDir/data/dgrp2dm6/reference/dmel-all-chromosome-r6.52.fasta"
        if not os.path.exists((reference + ".fai")):
            faidx = subprocess.run([samtools, "faidx", reference],
                                   capture_output=True, check=False)
    elif species == Species.HSAP.name:
        chrm_list = ["Nothing in here right now"]
        reference = "Add_humans_reference_latter"
    else:
        raise ValueError("Unknown species! Only Dmel and Hsap are supported!")

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

    # Prompt message
    # It might be good to have a progress bar in here
    print("Start processing the annotate VCF file. It might take a while...")

    with open(inputfile, "r", encoding="utf-8") as fp:
        current_block = []
        lines_in_block = []
        list_block = []
        list_lines_in_block = []
        previous_position = 0
        line_trck = 0
        # Read each line in the file and save to a list
        lines = []
        for line in fp:
            if line.startswith("##") or line.startswith("#"):
                continue

            # Extract variant information from the VCF fields
            fields = line.strip().split("\t")

            if len(fields) < 8:
                continue

            # Unpack FIELDS items
            # anc/ref, derived/alt
            chrom_str, pos, vcfid, ref, alt, qual, filtr, info = fields[
                :8
            ]  # maybe use pop() here

            # Cast pos to integer
            pos = int(pos)

            # Remove 'chr' from the chrom
            _, chrom = chrom_str.strip().split("r")

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

                refflank = flkng_bf + ref + flkng_af
                altflank = flkng_bf + alt + flkng_af

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
                if effect == "SYNONYMOUS_CODING" or effect == "NON_SYNONYMOUS_CODING":
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

                    if effect == "NON_SYNONYMOUS_CODING" or "SYNONYMOUS_CODING":
                        if siftscore == "NA":
                            deleteriousness = "NA"
                        else:
                            if float(siftscore) < 0.05:
                                deleteriousness = "deleterious"
                            else:
                                deleteriousness = "tolerated"
                    else:
                        deleteriousness = "NA"

            # Create a re-usable list to store the needed information
            line_list = [
                chrom_str,
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

    # Prompt message
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
        "-s",
        help="Species (either Dmel or Hsap)",
        dest="species",
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
    return parser


def main(argv):
    """
    This is the main program definition.
    """
    parser = parseargs()
    if argv[-1] == "":
        argv = argv[0:-1]
    args = parser.parse_args(argv)

    file = args.file
    outfile = args.outfile
    species = args.species
    nflankinbps = args.nflankinbps

    result = vcf_to_tsv(
        file=file, outfile=outfile, species=species, nflankinbps=nflankinbps
    )
    print(result)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        main(["-h"])
    else:
        main(sys.argv[1:])


# def main(argv):

#     starttime = time.time()
#     parser = parseargs()
#     if argv[-1] =='':
#         argv = argv[0:-1]
#     args = parser.parse_args(argv)


#     file = args.infile
#     outfile = args.outfile
#     species = args.species
#     nflankinbps = args.nflankinbps

#     # CMD for debugging
#     # file = "dgrp2dm6_chrX_rooted.vcf"
#     # outdirpath = None
#     # species = "Dmel"  #"reference/dmel-all-chromosome-r6.52.fasta"
#     # nflankinbps = 4

#     # Figure out where the input VCF files is
#     # and change the local path to it
#     path = os.getcwd()
#     workdir_, filename =  os.path.split(file)

#     # Here, it moves to the working directory
#     workdir = os.path.join(path, workdir_)
#     os.chdir(workdir)

#     # Define the basename for the outputs from the filename
#     basename, _ = filename.strip().split('.vcf')

#     # Include here steps for annotation with SNPeff and SIFT as a subprocess
#     # create the folders to save the output of each annotation step

#     # for SNPeff
#     snpeff_out = "annotations/snpeff"
#     if not os.path.exists(snpeff_out):
#               os.makedirs(snpeff_out)

#     # for SIFT
#     sift_out = "annotations/sift4g"
#     if not os.path.exists(sift_out):
#               os.makedirs(sift_out)

#     # Run SNPeff annotation
#     output_snpeff = (snpeff_out + "/" + basename + "_ann.vcf")
#     output_basename = (snpeff_out + "/" + basename)

#     # This control flows here allow the user to choose to bypass the SNP annotations if already done.
#     if doEFF:
#         SNPEFF="/Users/tur92196/snpEff/snpEff.jar"
#         SNPEFF_config="/Users/tur92196/snpEff/snpEff.config"
#         os.system("java -jar " + SNPEFF + " ann -c " + SNPEFF_config + " -fi " + "intervals/si_intervals.bed" + " -o gatk -csvStats " + output_basename +
#                   " -htmlStats " + (output_basename + ".html") + " -v " + dbSNPEFF + " " + (basename + ".vcf") + " > " + output_snpeff)

#     # Run SIFT annotation
#     if doSIFT:
#         SIFT4G = "/Users/tur92196/local/sift4g/SIFT4G_Annotator.jar"
#         os.system("java -jar " + SIFT4G + " -c -i " + output_snpeff + " -d " + dbSIFT + " -r " + sift_out)

#     if doTable:
#         # If None in the CMD was defined, used the same dir as input, otherwise check/use the one espeficied in 'outdir'
#         if outdir is None:
#             outdir = sift_out
#         else:
#             if not os.path.exists(outdir):
#                 os.makedirs(outdir)

#         outfile = (os.path.join(outdir, basename) + "_annotated_table.tsv")
#         header = ["chrom", "pos", "id", "ref", "alt", "refcount", "altcount","refflank", "altflank",
#                 "refcodon", "altcodon", "refaa", "altaa",
#                 "snpeff_effimpact", "snpeff_funclass", "snpeff_genename", "snpeff_trnscbiotype", "snpeff_genecoding", "snpeff_trnscid",
#                 "sift_trnscid", "sift_geneid", "sift_genename", "sift_region", "sift_vartype", "sifts_core", "sift_median", "sift_pred", "deleteriousness"]
#         header = ("\t".join(str(item) for item in header))

#         # This defines the output filename
#         output_sift = os.path.join(sift_out, (basename + "_ann_SIFTpredictions.vcf"))

#         # Prompt message
#         print("Start processing the annotate VCF file. It might take a while...")

#         with open(output_sift, 'r') as fp:
#             current_block = []
#             lines_in_block = []
#             list_block = []
#             list_lines_in_block = []
#             previous_position = 0
#             line_trck = 0
#             # Read each line in the file and save to a list
#             lines = []
#             for line in fp:
#                 if (line.startswith("##") or line.startswith("#")):
#                         continue

#                 # Extract variant information from the VCF fields
#                 fields = line.strip().split('\t')

#                 if len(fields) < 8:
#                         continue

#                 # Unpack FIELDS items
#                 # anc/ref, derived/alt
#                 chrom, pos, id, ref, alt, qual, filt, info = fields[:8] # maybe use pop() here

#                 # Cast pos to integer
#                 pos = int(pos)

#                 # Remove 'chr' from the chrom
#                 _, chrom = chrom.strip().split('r')

#                 # Unpack items in info
#                 altcount_, refcount_, snpeff_, *sift_ = info.strip().split(';')

#                 _, altcount = altcount_.strip().split('=')
#                 _, refcount = refcount_.strip().split('=')

#                 # Get the mutational context
#                 # Run samtools faidx as a python subprocess
#                 chrm_list = ["2L", "2R", "3L", "3R", "4", "X", "Y"]
#                 refflank, altflank = "NA", "NA"
#                 refcodon, altcodon, refaa, altaa = "NA", "NA", "NA", "NA"
#                 effect_impact, functional_class, gene_name, transcript_biotype, gene_coding, transcript_id = "NA", "NA", "NA", "NA", "NA", "NA"
#                 transcript, geneid, genename, region, varianttype, siftscore, siftmedian, siftpred, deleteriousness = "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"

#                 # Frim this part until the end, run only in canonical chromosomes
#                 if any(x == chrom for x in chrm_list):
#                     if not os.path.exists((reference + ".fai")):
#                         faidx = subprocess.run([samtools, "faidx", reference], capture_output=True)

#                     nt_bf = subprocess.run([samtools, "faidx", reference, (chrom + ":" + str(pos-nflankinbps) + "-" + str(pos-1) )], capture_output=True)
#                     nt_af = subprocess.run([samtools, "faidx", reference, (chrom + ":" + str(pos+1) + "-" + str(pos+nflankinbps) )], capture_output=True)

#                     _, flkng_bf, *_ = str(nt_bf.stdout).strip().split('\\n')
#                     _, flkng_af, *_ = str(nt_af.stdout).strip().split('\\n')

#                     refflank = (flkng_bf + ref + flkng_af)
#                     altflank = (flkng_bf + alt + flkng_af)

#                     # This part is for fixin basis in flaking bases string
#                     # if SNPs are close (less than nflankinbps)
#                     # This keeps track of block of positions with near SNPs
#                     current_position = pos

#                     if (current_position - previous_position) < nflankinbps: #(maybe nflankinbps - 1)
#                         if not current_block:
#                             current_block.append(previous_position)
#                             lines_in_block.append(line_trck-1)
#                         current_block.append(current_position)
#                         lines_in_block.append(line_trck)
#                     else:
#                         current_block = []
#                         lines_in_block = []

#                     if len(current_block) > 1:
#                         if current_block not in list_block:
#                             list_block.append(current_block)
#                         if lines_in_block not in list_lines_in_block:
#                             list_lines_in_block.append(lines_in_block)

#                     previous_position = current_position

#                     # Unpack items in SNPeff info::eff
#                     _, eff_ = snpeff_.strip().split('=')
#                     effect, eff  = eff_.strip().split('(')

#                     # ALL SNPEff fields
#                     effect_impact, functional_class, codon_change, aa_change, gene_name, transcript_biotype, gene_coding, transcript_id, *exon  = eff.strip().split('|')

#                     # Unpack codons
#                     if (effect == "SYNONYMOUS_CODING" or effect == "NON_SYNONYMOUS_CODING"):
#                         refcodon, altcodon = codon_change.strip().split('/')

#                     # Unpack items in info::sift (when present)
#                     if len(sift_) > 0:
#                         # Unpack items in info::eff
#                         _, siftinfo = sift_[0].strip().split('=')

#                         # ALL SIFT fields
#                         allele, transcript, geneid, genename, region, varianttype, aa, aaposition, siftscore, siftmedian, siftnumseqs, alleletype, siftpred =   siftinfo.strip().split('|')
#                         refaa, altaa = aa.strip().split("/")

#                         if effect == "NON_SYNONYMOUS_CODING" or "SYNONYMOUS_CODING":
#                             if  siftscore == 'NA':
#                                 deleteriousness = 'NA'
#                             else:
#                                 if float(siftscore) < 0.05:
#                                         deleteriousness = "deleterious"
#                                 else:
#                                     deleteriousness = "tolerated"
#                         else:
#                             deleteriousness = 'NA'

#                 # Create a re-usable list to store the needed information
#                 line_list = [chrom, pos, id, ref, alt, refcount, altcount, refflank, altflank, refcodon, altcodon, refaa, altaa, effect_impact, functional_class, gene_name, transcript_biotype, gene_coding, transcript_id, transcript, geneid, genename, region, varianttype, siftscore, siftmedian, siftpred, deleteriousness]
#                 lines.append(line_list)
#                 line_trck +=1

#         # Prompt message
#         # print("Breaking larger blocks of SNPs...")
#         # new_list_block = []
#         # new_list_lines_in_block = []

#         # for block_lines, block_pos in zip(list_lines_in_block, list_block):
#         #     #print(block_pos)
#         #     sub_block = [block_pos[0]]
#         #     sub_block_lines = [block_lines[0]]

#         #     for i in range(1, len(block_pos)):
#         #         distance = block_pos[i] - sub_block[0]  # Calculate distance from the first element

#         #         if distance <= 3:
#         #             sub_block.append(block_pos[i])
#         #             sub_block_lines.append(block_lines[i])
#         #         else:
#         #             new_list_block.append(sub_block)
#         #             sub_block = [block_pos[i]]

#         #             new_list_lines_in_block.append(sub_block_lines)
#         #             sub_block_lines = [block_lines[i]]

#         #     # Add the last sub-block
#         #     new_list_block.append(sub_block)
#         #     new_list_lines_in_block.append(sub_block_lines)

#         # Prompt message
#         print("Fixing ancestral and derived mutations on flanking bases. It might take a while...")
#         #for block_lines, block_pos in zip(new_list_lines_in_block, new_list_block):
#         for block_lines, block_pos in zip(list_lines_in_block, list_block):
#             blocKsize = len(block_lines)
#             if blocKsize > 1:
#                 for i, current_element in enumerate(block_lines):
#                     elements_before = block_lines[:i]
#                     elements_after = block_lines[i + 1:]
#                     current_pos = block_pos[i]
#                     pos_before = block_pos[:i]
#                     pos_after = block_pos[i + 1:]

#                     # Get the flanking bases that need to be updated
#                     current_refflanking = list(lines[current_element][7])
#                     current_altflanking = list(lines[current_element][8])

#                     # Check if there are multiple elements before or after
#                     if len(elements_before) >= 1:
#                         for j, eb in enumerate(elements_before):
#                             dist = abs(current_pos - pos_before[j])
#                             if dist > nflankinbps:
#                                 continue
#                             else:
#                                 ref_before = (lines[eb][3])
#                                 alt_before = (lines[eb][4])
#                                 current_refflanking[nflankinbps-dist] = ref_before
#                                 current_altflanking[nflankinbps-dist] = alt_before

#                     if len(elements_after) >= 1:
#                         for j, ea in enumerate(elements_after):
#                             dist = abs(current_pos - pos_after[j])
#                             print(dist)
#                             if dist > nflankinbps:
#                                 continue
#                             else:
#                                 ref_after = (lines[ea][3])
#                                 alt_after = (lines[ea][4])
#                                 current_refflanking[nflankinbps+dist] = ref_after
#                                 current_altflanking[nflankinbps+dist] = alt_after

#                     # Here update the flanking sequences in the corresponding line
#                     lines[current_element][7] = "".join(current_refflanking)
#                     lines[current_element][8] = "".join(current_altflanking)

#         # Prompt message
#         print("Exporting the processed VCF to a .TSV file")
#         with open(outfile, 'a') as fo:
#             fo.write(header + "\n")
#             for line in lines:
#                 fo.write(("\t".join(str(item) for item in line)) + "\n")
#         print("DONE!!!")


# if __name__ == "__main__":

#     if len(sys.argv) < 2:
#         main(['-h'])
#     else:
#         main(sys.argv[1:])
