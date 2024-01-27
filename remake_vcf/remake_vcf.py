'''
Remake a vcf to a standard format
'''

import subprocess  # this is being used to run samtools
from pathlib import Path  # This is being used to check if the file exists

# This is experimental
# class VcfRemaker:
#     def __init__(self, input_vcf_file: str, output_vcf_file: str, chrom_name: str, path_to_reference: str = "reference") -> None:
#         self.input_vcf_file = input_vcf_file
#         self.chrom_name = chrom_name

#     @classmethod
#     def run_faidx(cls, path_to_reference: str, chrom_name: str) -> None:
#         samtools_exec = "/usr/local/anaconda3/envs/bioinfo/bin/samtools"
#         subprocess.run([samtools_exec, "faidx", f"{path_to_reference}/dmel-{chrom_name}-chromosome-r5.57.fasta"], capture_output=True, check=True)

# @classmethod
# def remake_vcf(cls, vcf_file: str, create_faidx: bool, path_to_reference: str = None, rename_chrom: bool, chrom_name: str, output_file: str):
#     if create_faidx:
#         if path_to_reference is None:
#             raise ValueError("path_to_reference must be provided if create_faidx is True")
#         elif chrom_name is None:
#             raise ValueError("chrom_name must be provided if create_faidx is True")
#         else:
#             cls.run_faidx()  # Call the run_faidx method

#     if rename_chrom:
#         if chrom_name == "":
#             raise ValueError("chrom_name must be provided if rename_chrom is True")
#     else:
#         raise ValueError("chrom_name must not be provided if rename_chrom is False")


# def generate_file_names(population_prefix: str, chromosome_name: str, input_file_path: str = "", reference_file_path: str = "reference") -> Tuple[str, str]:
#     """
#     Generate input and reference file names
#     """
#     input_file = Path(input_file_path) / f"{population_prefix}_chr{chromosome_name}.vcf"
#     reference_file = Path(reference_file_path) / f"dmel-{chromosome_name}-chromosome-r5.57.fasta"
#     return input_file, reference_file


# create faidx with samtools if not exists
def create_samtools_faix(reference_file: Path, samtools_path: Path) -> None:
    """
    Check the existence of a reference file faidx
    or create it with samtools
    """

    samtools = samtools_path
    reference_fai_file = reference_file.with_suffix(reference_file.suffix + ".fai")
    if not reference_fai_file.exists():
        make_faidx = subprocess.run([samtools, "faidx", reference_file], capture_output=True, check=True)


def remake_vcf_header(input_file: Path, reference_file: Path, chrom_name: str, output_file: Path) -> None:
    """
    Remake a vcf header to a standard format compatible with GATK
    """
    # Process header lines
    header_lines = []
    with open(input_file, 'r', encoding='utf-8') as input_file_obj:
        # Process lines starting with "##" or "#"
        for line in input_file_obj:
            if line.startswith("##") or line.startswith("#"):
                header_lines.append(line)

    # Add new lines and modify existing lines
    modified_header_lines = header_lines.copy()  # Create a copy of the original header lines

    # Remove two lines but add then later in the file processing
    chrm_line_info = modified_header_lines.pop(1)
    format_str_inf = modified_header_lines.pop(1)
    colum_names_line = modified_header_lines.pop(1)

    # Add new line
    modified_header_lines.append("##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">\n")

    # Add FORMAT line remove above
    modified_header_lines.append(format_str_inf)

    # Add multiple new lines
    new_lines = [
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n"
    ]
    modified_header_lines.extend(new_lines)

    # Add a new and modified contig/chromosome line information
    modified_header_lines.append(f"##contig=<ID={chrom_name},length=23011544>\n")

    # Add line with reference value
    modified_header_lines.append(f"##reference=file://{reference_file}\n")

    # Add column names line
    modified_header_lines.append(colum_names_line)

    # Write the modified header lines to the output file or perform further manipulations
    # For example, writing to the output file:
    with open(output_file, 'w', encoding='utf-8') as output_file_obj:
        for line in modified_header_lines:
            output_file_obj.write(line)



def remake_vcf(input_file: Path, reference_file: Path, chrom_name: str, output_file: Path, samtools_path: Path) -> None:
    """
    Remake a vcf to a standard format compatible with GATK
    """
    
    # Check the existence of a reference file faidx
    create_samtools_faix(reference_file, samtools_path)
    print("Faidx created...")

    # Run remake_vcf_header
    remake_vcf_header(input_file, reference_file, chrom_name, output_file)
    print("New header created...")

    # Process remaining lines
    with open(input_file, "r", encoding="utf-8") as input_file_obj, open(output_file, "a", encoding="utf-8") as output_file_obj:
        for line in input_file_obj:
            if (line.startswith("##") or line.startswith("#")):
                continue

            # Extract variant information from the VCF fields
            fields = line.strip().split("\t")

            if len(fields) < 9:
                continue

            # Fix chromosome name
            if fields[0].startswith("chr"):
                chrom_str = fields[0]
            elif fields[0] == "." or fields[0] == "1":
                fields[0] = "chr" + chrom_name
                chrom_str = fields[0]
            else:
                fields[0] = "chr" + fields[0]
                chrom_str = fields[0]

            # Get vcf ref and alt bases
            vcf_candidate_ref_base = fields[3]
            vcf_candidate_alt_bases = fields[4].split(",")

            # Create a list with all bases: reference and alternatives
            # The assigned snps-sites reference should the the first
            vcf_candidate_bases = [vcf_candidate_ref_base] + vcf_candidate_alt_bases

            # Create a list of indexes based on the original order of vcf_candidate_bases elements
            vcf_candidate_bases_idx = [i for i, x in enumerate(vcf_candidate_bases)]

            # Check reference base with samtools faidx
            pos = int(fields[1])  # First, cast pos to integer
            _, chrom = chrom_str.strip().split("r")  # Remove 'chr' from the chrom

            # Call samtools faidx to check reference base
            genome_ref_base = subprocess.run([
                samtools_path,
                "faidx",
                reference_file,
                (chrom + ":" + str(pos) + "-" + str(pos)),
            ],
                                    capture_output=True,
                                    check=False)

            _, genome_ref_base, *_ = str(genome_ref_base.stdout).strip().split("\\n")

            # Get the index of the missing data character '*'
            idx_of_missing_data = vcf_candidate_bases.index("*")

            # Convert the missing data code to the standard format
            fields[9:] = ["./." if x == str(idx_of_missing_data) else x for x in fields[9:]]

            # Then remove the missing data character index from vcf_candidate_bases_idx
            vcf_candidate_bases_idx.remove(idx_of_missing_data)

            # Get the index of the genome reference base on the list of vcf candidate bases
            idx_of_ref_base = vcf_candidate_bases.index(genome_ref_base)

            # Create a list of re-ordered indexes where the first value is the index of the reference base
            vcf_candidate_bases_idx.sort(key=lambda x: x != idx_of_ref_base)

            # Convert the genotypes to standard format
            # and count the number of each non-missing genotype
            genotype_counts = []
            for idx, idxval in enumerate(vcf_candidate_bases_idx):
                print(idx, idxval)
                fields[9:] = [f"{idx}/{idx}" if x == str(idxval) else x for x in fields[9:]]
                genotype_counts.append(fields[9:].count(f'{idx}/{idx}'))

            # New list of vcf bases after re-ordered and reference base identified
            vcf_final_bases = [vcf_candidate_bases[i] for i in vcf_candidate_bases_idx]

            # Replace the old vcf reference and alternative bases with the new ones
            fields[3] = vcf_final_bases[0]
            vcf_final_bases.pop(0)

            fields[4] = ",".join(vcf_final_bases)

            # Save reference counts and frequencies
            total_genotype_count = sum(genotype_counts)
            ref_counts_value = genotype_counts.pop(0)
            ref_counts = f"RC={ref_counts_value}"
            ref_frequencies = f"RF={ref_counts_value/total_genotype_count}"
            
            # Save alternative counts and frequencies
            alt_counts = f"AC={','.join(str(x) for x in genotype_counts)}"
            alt_frequencies = f"AF={','.join(str(x / sum(genotype_counts)) for x in genotype_counts)}"
            fields[7] = alt_counts + ";" + alt_frequencies   

            # Write the new line
            reassembled_line = '\t'.join(fields)

            with open(output_file, 'a', encoding='utf-8') as output_file_obj:
                output_file_obj.write(reassembled_line + "\n")


def main() -> None:

    # Define input and output files
    input_file = "ZI_test_chr2L.vcf"
    reference_file = Path("reference/dmel-2L-chromosome-r5.57.fasta")
    chrom_name = "2L"
    output_file = "ZI_test_chr2L_reassembled.vcf"
    samtools_path = Path("/usr/local/anaconda3/envs/bioinfo/bin/samtools")
    remake_vcf(input_file, reference_file, chrom_name, output_file, samtools_path)


if __name__ == "__main__":
    main()