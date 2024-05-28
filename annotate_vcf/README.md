# Annotation pipeline

Here I show how to run the annotation pipeline.
```zsh
python pipeline_to_annotate_vcf.py -i examples/lifted/example_remade_rooted_lifted.vcf -d Drosophila_melanogaster -b examples/intervals/dm6_short_introns.bed -o examples/snpeff -s /Users/tur92196/local/sift4g/BDGP6.83 -f examples/sift4g
```

To run this pipeline like above, you should have installed [SNPEff](https://pcingola.github.io/SnpEff/) and [SIFT4](https://sift.bii.a-star.edu.sg/sift4g/). SIFT4 is optional if you don't provide a PATH to a database. But if you provided, you should provid the PATH for an output folder. Another optional argument is the PATH for the file containing intervals you want to include in SNPEff annotation. It should be a BED-like file containing somehow "custom" annotations. The pipeline looks for the presence of a file PATH and when triggered, it implements SNPEff `-interval` argument.

The other parameters are self explained when you have the pipeline help message:
```zsh
usage: python pipeline_to_annotate_vcf.py [-h] -i INPUT_FILE -d SNPEFF_DATABASE [-b SNPEFF_INTERVAL_FILE] -o SNPEFF_OUTPUT_FOLDER [-s SIFT4G_DATABASE]
                                          [-f SIFT4G_OUTPUT_FOLDER]

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE         input vcf file not annotated (default: None)
  -d SNPEFF_DATABASE    SNPEff database (default: None)
  -b SNPEFF_INTERVAL_FILE
                        SNPEff interval BED file (default: None)
  -o SNPEFF_OUTPUT_FOLDER
                        output folder for SNPEff annotations (default: None)
  -s SIFT4G_DATABASE    SIFT4G database (default: None)
  -f SIFT4G_OUTPUT_FOLDER
                        Output folder for SIFT4G annotations (default: None)
```

### Issues:
- `simplify_snpeff.py` only takes the first term from SNPEff annotation. For bi-allelic SNPs, the first term should be the important one (with the annotation of the most impacted change); but for tri-allelic there are the issue with the most important annotation per allele, and which allele is the most (or least) common one. Right now, the script picks the first term, regardless of all these issues for tri-allelic(s). ADVICE: remove tri-allelics before running this pipeline.
- `vcf_to_tsv.py` also had to deal with tri-allelic SNPs in SIFT4g annotations. I did a quick fix to baypass it. It basica disregard the second annotated term, basically discarding the functional annotation of the second alternative allele. Because of it, the logic would be to remove tri-allelic SNPs.
- `get_reversed_complementary_strand()` also doesn't know how to deal with tri-allelics.