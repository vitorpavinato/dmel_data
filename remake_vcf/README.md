# Remake vcf file

Some ideas to make DGN data usable for vcf-based tools.

I recommend you structure your working directory as something like this:
```bash
.
├── dpgp3
│   ├── dmel_data
│   ├── dpgp3.txt
│   ├── masked
│   │   ├── Chr2L
│   │   ├── Chr2R
│   │   ├── Chr3L
│   │   ├── Chr3R
│   │   ├── fastas
│   │   └── vcfs
│   │       ├── liftedover
│   │       ├── remade
│   │       └── rooted
│   └── originals
├── reference
├── packages
└── simulans_sequences
```

### Download DGN data
First download the consensus sequence FASTA from [DGN](https://www.johnpool.net/genomes.html) for the target population (sample). Go to the bottom of the page where you see a bunch of links like \<name\>_SEQ. These are the links to a compressed folder that contains `.seq` files with consensus FASTA sequences for each individual organized for each chromosome. Here I am showing as an example the analysis of data from DPGP3 that contains the 197 Zambia genomes.

```zsh
mkdir data
cd data
wget http://pooldata.genetics.wisc.edu/dpgp3_sequences.tar.bz2

md5sum-lite dpgp3_sequences.tar.bz2 906d282740a56e5273a4cbc5abfb61f9
```

After checking the md5sum and decompressing the file, you should place the folder containing the chromosomes FASTA inside the folder I named `originals`.

### Run DGN masking scripts
Second, download the [masking package](http://johnpool.net/masking.zip) also from DGN. It contains a set of files and scripts used to mask problematic sites previously identified containing traces of identity-by-descent (IBD) or admixture. Copy both scripts and the interval files (two `.csv` files) inside each chromosome folder (that should now be inside `originals`), and run the two `Perl`` scripts. Here I am showing how to run them inside chromosome 2L folder (I am showing only the analysis of Chr2L, for the rest of this document).

First download the masking package:
```zsh
wget http://johnpool.net/masking.zip
```

You can either delete the compressed files or keep them. You can create a folder in the root of your project called `packages` to save the scrips, CSVs and the original `.zip` files in there.

Then copy these files to each of the folder containing sequences for each chromosome.
```zsh
cd chr2L

# Run the IBD masking script first:
perl ibd_mask_seq.pl

# Wait until it finish, then run the Admixture masking script:
perl admixture_mask_seq.pl
```

You need to repeat it for all folders (I know, it is tedious, I stopped at chr3R after running the other three major chromosomes - maybe you can open four terminal windows and run then one script at time, but in four different chromosomes...)

If the individual has IBD or Admixture traces, an `N` will replace the nucleotide on the problematic position. For each chromosome, the result is a set of individuals genomes (masked or not). Move the masked files to the corresponding folder named with the chromosome name inside the folder `masked`.

```zsh
mv *_Chr2L.seq ../../masked/Chr2L/
```

### Add sample names to each sequence
These `.seq` files don't have header (I don't know why?). Copy the name of each individual and place it as a header (using bash please with something like ```ls -1 *.seq | rev | cut -c11- | rev > dpgp3.txt```). You can keep this files at the root directory of your project, as it is shown in the above directory structure. I prepared a script that takes the names of each individual placed in a one-column `.txt` file (as the one created above), opens the corresponding `.seq`, and adds the corresponding name of the sample as the header:

```zsh
# You should run it for each chromosome file:
cd ../..masked/Chr2L
python ../../dmel_data/remake_vcf/add_head_to_fasta_files.py -s _Chr2L.seq -n ../../dpgp3.txt
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

A note here: If you are trying to install snp-site from conda on a Apple computer with Mx chips, run this command before (with the target conda environemtn activated):
```zsh
conda config --env --set subdir osx-64
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
JAVA='/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home/bin/java'
GATK='/Users/tur92196/local/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar'
```

To run GATK `LiftoverVcf`, you should have:
- your `.vcf` you want to lift to the newest genome;
- the reference sequence of the target genome (the one you are lifting the SNPs to);
- the associated `.dict` file (see how to get one below);
- and the chain files, downloaded from UCSC (found [here](https://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/)), containing the chains from the oldest (here dm3) to the newest target genome (dm6).

Go inside the folder you have the genome file of the target genome. There we will run, to obtain the `.dict` file. You can find the [documentation here](https://gatk.broadinstitute.org/hc/en-us/articles/360037422891-CreateSequenceDictionary-Picard-):

```bash
cd reference
$JAVA -jar $GATK CreateSequenceDictionary \\ 
    -R dm6.fa \\ 
    -O dm6.fa.dict
```

With everything in place, then type, to have lift SNPs from one genome to the target genome:
```bash
cd ..

$JAVA -jar $GATK LiftoverVcf \\ 
    -I dmel_data/remake_vcf/example/example_output_remade_rooted.vcf \\ 
    -O dmel_data/remake_vcf/example/example_output_remade_rooted_lifted.vcf \\ 
    -C reference/dm3ToDm6.over.chain \\ 
    --REJECT dmel_data/remake_vcf/example/example_output_remade_rooted_rejected.vcf \\ 
    -R reference/dm6.fa \\ 
    --WARN_ON_MISSING_CONTIG true \\ 
    --RECOVER_SWAPPED_REF_ALT true &
```

### Additional filtering

Now that everthing is in place, you can filter out SNPs with more than one alternative allele (retain only bi-allelic), and SNPs with a lot of missing data (more than 50%). We are going to use `vcftools` for this task.

```zsh
vcftools --vcf remake_vcf/example/example_output_remade_rooted_lifted.vcf \\ 
        --out example_output_remade_rooted_lifted_fltr  
        --min-alleles 2 \\ 
        --max-alleles 2 \\ 
        --max-missing 0.5 \\ # Depending what you want to accomplish, this might be too restrictive
        --recode \\ 
        --recode-INFO-all
```

Remove the `.recode.` part of the file name to kept file naming simple.
```zsh
mv example_output_remade_rooted_lifted_fltr.recode.vcf example_output_remade_rooted_lifted_fltr.vcf
```