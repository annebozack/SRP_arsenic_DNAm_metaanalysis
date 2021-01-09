# Exposure to arsenic at different life-stages and DNA methylation meta-analysis in buccal cells and leukocytes 

This repository contains the necessary scripts to reproduce the analysis of
Bozack et al.'s "Exposure to arsenic at different life-stages and DNA
methylation meta-analysis in buccal cells and leukocytes". A preprint of the
manuscript can be found here **(insert link)**.

The organization of the repository is as follows:

- The `bangladesh-study` folder contains the scripts and results associated
  with the Bangladesh DNAm studies.
- The `chile-study` directory is made up of subdirectories containing the code
  and results associated with the buccal cell and the PBMC DNAm analyses. It
  also holds a  directory with notebooks assessing the within-participant
  buccal-PBMC sample similarities, and a directory with a notebook of
  descriptive statistics used to create Table 1 in the accompanying paper.
- The `DMR-meta-analysis` folder contains the `DMR-meta-analysis.Rmd` notebook,
  which details the meta-analysis procedure and summarized the results output
  by [comb-p](https://github.com/brentp/combined-pvalues).
- The `helper-scripts` directory contains multiple helper files used to perform
  the analyses described in the paper.
