# How to obtain SFSs from tsv files

I made available in the `notebooks` folder two notebooks I created whan I was analyzing *Drosophila melanogster* data. It was easy to report the steps of getting the SFSs from the data while developing some code support for the analysis. This part of the repository is the most tested with a lot of scripts implementing pytest(s). Functions were organized in which context they might apply, for example `sfs_utils.py` contain all code related to obtaining SFS(s). These notebooks are here for completeness and reference. This repository should not have any updated especially because now it is associated with the manuscript:

```dotnetcli
"Poisson random field ratios in population genetics: estimating the strength of selection while sidestepping non-selective factors" by Jody Hey and Vitor Pavinato.
```

Here is the link to the repository containing the program we developed called [PRF_ratios](https://github.com/vitorpavinato/PRF_Ratios).

### To do's:
- [ ] Make pytest for snp_utils.py:
    - [ ] filter_short_introns_from_bed()
    - [x] filter_snps_by_interval()
    - [x] create_snp_total_counts_dict()
