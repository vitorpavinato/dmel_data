"""
This module contains the function to fix mutational context
of processed vcf lines.
"""

from typing import List
from process_vcf_lines_utils import get_reversed_complementary_strand


def fix_mutational_context(
        list_lines_in_block: List[List[int]], 
        list_block: List[int], 
        vcf_lines: List[List[str]], 
        nflankinbps: int
) -> List[List[str]]:
    """
    This function fix the mutational context of processed vcf lines.
    It runs at the end, when all lines were processed. It fixes the
    refcontext and altcontext of closest SNPs.
    """

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

    return vcf_lines