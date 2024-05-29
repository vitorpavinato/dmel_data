## Drosophila melanogaster data (Dmel_data)

This repository documents the steps to obtain a pair of Site-frequency Spectrums to test our method PRFratio for the estimation of $2Ns$.

Table of Content:

- [Remake vcf file](#Remake-vcf-file)
- [Annotate vcf SNPs](#Annotate-vcf-SNPs)
- [Create SFSs sets](#Create-SFSs-sets)


### Remake vcf file
Follow the steps on the `remake_vcf` page to make DGN files usable (it also have the steps to lift dm3 SNPs to dm6 coordinates).

### Annotate vcf SNPs
A pipeline and jupiter notebooks in `annotate_vcf` folder that documents all the steps to annotate `.vcf` files.

### Create SFSs sets
Jupyter notebooks in `create_sfss` folder with all steps to get the pairs of SFSs.
