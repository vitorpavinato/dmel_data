# Better VCF from Drosophila Genome Nexus (DGN)

Some ideas to make DGN data usable for vcf-based tools.

### Download DGN data
First download the consensus sequence FASTA from [DGN](https://www.johnpool.net/genomes.html). Go to the bottom of the page where you see a bunch of links like \<name\>_SEQ. These are the links to a compressed folder that contains `.seq` files with consensus FASTA sequences for each individual organized for each chromosome. We used the data from DPGP3 that contains the 197 Zambia genomes.

### Run DGN masking scripts
Second, download the [masking package](http://johnpool.net/masking.zip) also from DGN. It contains a set of files and scripts used to mask problematic sites previously identified containing traces of identity-by-descent (IBD) or admixture. Copy both scripts and the interval files (two `.csv` files) inside each chromosome folder, and run the two `Perl`` scripts. Here I am showing how to run them inside chromosome 2L folder:

```zsh
cd chr2L

# Run the IBD masking script first:
perl ibd_mask_seq.pl

# Wait until it finish, then run the Admixture masking script:
perl admixture_mask_seq.pl
```

You need to repeat it for all folders (I know, it is tedious, I stopped at chr3R after running the other three major chromosomes - maybe you can open four terminal windows and run then one script at time, but in four different chromosomes...)

If the individual has IBD or Admixture traces, an `N` will replace the nucleotide on the problematic position. For each chromosome, the result is a set of individuals genomes (masked or not). 

### Add sample names to each sequence
These `.seq` files don't have header (I don't know why?). Copy the name of each individual and place it as a header (using bash please). I prepared a script that takes the names of each individual placed in a one-column `.txt` file, opens the corresponding `.seq`, and add the corresponding name of the sample as the header:

```zsh
python add_head_to_fasta_files.py -s _Chr2L.seq -n dpgp3.txt
```

### Run snp-site
Then concatenate all `.seq` file by chromosome using `cat`. I did something like this. Inside the folder with all sequences from the same chromosome, type:

```zsh
cat *.seq > ZI_Chr2L.fasta
```

Then, run [snp-site](https://sanger-pathogens.github.io/snp-sites/) to convert the *multi-fasta alignment* to a somehow descent and usable `.vcf` file. Here is how the command looks like:

```zsh
snp-sites -v -o ZI_Chr2L.vcf ZI_Chr2L.fasta
```

This will create a `.vcf` from the `.fasta` file. Here is important to have all sequences with their corresponding header. Usually it is the sample name. The VCF will have then the same ID as the FASTA for each sample. 

Snp-site has no way to know which base in the multi-fasta alignment is in the reference base, so I prepared a script to deal with this (see [remake_vcf](https://github.com/vitorpavinato/dmel_data/tree/main/remake_vcf) folder). The script also identify which genotype number corresponds to the missing character `*` and change the number to `./.`, effectively making the masking done before usable. And it re-sort the order of alternative alleles and make genotype codes accordingly to the order. The most important ordering is the true reference and the alternatives, so the allele that is the reference in the genome (dm3 or the flybase release 5) will be `0` and reference genotypes will be `0/0`. 

Here is how you run `remake_vcf`.
```zsh
python remake_vcf.py -i example.vcf -r reference/dm3.fa -c chr2L -l 23011544 -o example_output.vcf -s /usr/local/anaconda3/envs/bioinfo/bin/samtools
```

You need to provid the path for samtools, the name of the chromosome and the length (length is only used to annotate the final `.vcf` file).

I advise to use the reference genomes (and annotations) stored at [UCSC genome browse](https://genome.ucsc.edu/cgi-bin/hgGateway) because the chromosomes are named with the same naming system as in the chain file. The chain files will be used in later steps to liftover the SNPs from dm3 (dmel r.5.x) to dm6 (dmel r.6.x).

### Liftover VCF
With a `better DGN vcf file`, now is time to lift the positions over to the newest (at least the newest at UCSC genome browser) Dmel genome (release 6). We are going to use GATK `Picard` [LiftoverVcf](https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard). 

GATK is a great tool, not easy to use. Sometimes the headaches that follows when you try to use this tools is caused by some incompatibilities with the Java you have installed. To avoid it, we are going to use GATK packaged in a container. First, you need to have Docker or Singularity installed (sorry, I am not going to cover how to install here). Then go to this [page](https://hub.docker.com/r/broadinstitute/gatk/tags/) to find the version of GATK Docker you want to have. With the version, and after installing Docker in your computer, open a terminal window and type:

```zsh
docker pull broadinstitute/gatk:latest
```
This will pull a copy of the container in your computer (it is like a virtual machine packaged with everything you need to run the target program).

Now, you need to be able to use the tools inside the Docker with your data. There are mainly two ways to do so: one is to copy your data to a folder inside the machine, the second and the one showed here, is to link a native folder in your computer to the folder inside the container. You machine folder, then should contain all the files you need for the liftover:
- your target `.vcf`;
- the reference sequence of the target genome (the one you are lifting the SNPs to);
- the associated `.dict` file (see how to get one below);
- and the chain files, downloaded from UCSC.

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

With everything in place, then type, to have lift SNPs from one genome to the target genome:
```bash
cd ..

gatk LiftoverVcf \  
    I=example/example_output.vcf \ 
    O=liftover/example_output_lifted.vcf \  
    CHAIN=liftover/dm3ToDm6.over.chain \  
    REJECT=example_output_rejected.vcf \  
    R=liftover/dm6.fa

```