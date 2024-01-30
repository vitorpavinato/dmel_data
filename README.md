## Better VCF from Drosophila Genome Nexus (DGN)

First download the consensus sequence fasta from [DGN](https://www.johnpool.net/genomes.html). We used the data from DPGP3 that contains the 197 Zambia genomes.

Second, download the [masking package](http://johnpool.net/masking.zip) also from DGN. It contains a set of files and scripts used to mask problematic sites with traces of identity-by-descent (IBD) and admixture. Run the two Perl scripts inside the folder containing the sequence data (each individual genome data is stored by chromosome). If the individual has IBD or Admixture traces, an N will replace the nucleotide in that position. For each chromosome, the result is a set of individuals genomes (masked or not). These fasta don't have header (I don't know why?). Copy the name of each individual and place it as a header (usin bash please). Then concatenate all sequences together by chromosome.

Then, run [snp-site](https://sanger-pathogens.github.io/snp-sites/) to convert the multi-fasta alingment to a (decent) vcf file. Snp-site has no way to know which base is in the reference, so I prepared a script to deal with it (see [remake_vcf](https://github.com/vitorpavinato/dmel_data/tree/main/remake_vcf) folder). The script also identify which genotype number corresponds to the missing character `*` and change the number to `./.`, effectively making the masking done before usable. And it re-sort the order of alternative alleles and make genotype codes accordly to the order. The most important ordering is the true reference and the alternatives, so the allele that is the reference in the genome (dm3 or the flybase release 5) will be `0` and reference genotypes will be `0/0`.

However, I advice to use the reference genomes (and annotations) stored at [UCSC genome browse](https://genome.ucsc.edu/cgi-bin/hgGateway) because the chromosomes are named with the same naming system as in the chain file. The chain files will be used in later steps to liftover from dm3 (dmel r.5.x) to dm6 (dmel r.6.x).

### Liftover VCF
With a `better DGN vcf file`, now is time to lift the positions over to the newest (at least the newest at UCSC genome browser) Dmel genome (release 6). We are going to use GATK Picard [LiftoverVcf](https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard). 

GATK is a great tool, not easy to use. Sometimes the headaches that follow when you try to use this tools is caused by some incompatibilities with the Java you have installed. To avoid it, we are going to use GATK packed in a container. First, you need to have Docker or Singularity installed (sorry, I am not going to cover how to install here). Then go to this [page](https://hub.docker.com/r/broadinstitute/gatk/tags/) to find the version of GATK Docker you want to have. With the version, and after installing Docker in your computer, open a terminal window and type:

```zsh
docker pull broadinstitute/gatk:latest
```
This will pull a copy of the container in your computer (it is like a virtual machine packaged with everything you need to run the target program).

Now, you need to be able to use the tools inside the Docker with your data. There mainly two ways to do so, one is to copy your data to a folder inside the machine, the second and the one showed here, is to link a native folder in your computer to the folder inside the container. This should contain all the files you need for the liftover:
- your target vcf;
- the reference sequence of the target genome (the one you are lifting the snps over);
- the associated .dict file (see how to get one below);
- and the chain files, downloaded at UCSC.

How to create a folder to link the machines:
```zsh
docker run -v ~/Documents/Repositories/dmel_data/remake_vcf:/gatk/my_data -it broadinstitute/gatk:latest
```

Where `my_folder` is a folder in my computer and `gatk/my_data` is a folder inside the Docker.

Put everything you have except the `.dict` file in your folder. After creating the link, we will be able to navigate inside the container. There we will run, to obtain the `.dict` file:

```bash
cd liftover
gatk CreateSequenceDictionary \ 
    R=dm6.fa \ 
    O=dm6.fa.dict
```

With everything in place, then type:
```bash
cd ..

gatk LiftoverVcf \  
    I=example/example_output.vcf \ 
    O=liftover/example_output_lifted.vcf \  
    CHAIN=liftover/dm3ToDm6.over.chain \  
    REJECT=example_output_rejected.vcf \  
    R=liftover/dm6.fa

```