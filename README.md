# dePCR Sequencing Workflow

This repository contains a reproducible workflow for processing dePCR sequencing data, generating primer-template count matrices, exploring primer-template composition, and fitting Bayesian distributional negative binomial models with covariate-dependent overdispersion.

## Repository description

The workflow starts with publicly available sequencing data from NCBI SRA and proceeds through paired-end read merging, quality control, primer-template feature counting, metadata integration, exploratory analysis, Bayesian model fitting, model comparison, and exploratory transfer assessment.

The main objective is to evaluate how primer sequence configuration, mismatch status, annealing temperature, and experimental design affect both expected primer-template abundance and residual count variability.

## Data sources

The workflow uses sequencing data from the following NCBI BioProjects:

- [PRJNA513137](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA513137)
- [PRJNA1072695](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1072695)

The analysis also uses:

- SRA run-information tables for both BioProjects;
- Supplementary sample metadata distributed with the original study;
- A primer-template feature annotation table;
- An experiment-specific feature table defining the expected features for each experiment;
- A feature metadata table containing mismatch status and nucleotide identity at the evaluated primer positions.

An example feature annotation file is provided at:

[`files/list_templates.txt`](files/list_templates.txt)

Expected columns:

```text
feature    template    primer
