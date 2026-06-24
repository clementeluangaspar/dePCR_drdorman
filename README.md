# Data thinning and Negative Binomial with covariate-dependent overdispersion for dePCR dataset.
# DePCR Sequencing Workflow

This repository contains a reproducible workflow for processing DePCR sequencing data, generating primer-template count matrices, applying data thinning, and fitting negative binomial models with covariate-dependent overdispersion.

## Repository description

Reproducible DePCR sequencing workflow for preprocessing, quality control, primer-template counting, data thinning, negative binomial modeling, and generation of figures and tables.

## Overview

I use this workflow to process sequencing data from the DePCR experiment. The analysis starts with raw SRA files, converts them to FASTQ, merges paired-end reads, performs quality control, converts merged reads to FASTA/FNA format, and counts primer-template combinations.

After preprocessing, I use the count matrix to evaluate sequencing and merging metrics, apply data thinning, and model counts using negative binomial regression with covariate-dependent overdispersion.

## Workflow organization

I divide the analysis into five RMarkdown files. Each file documents one major step of the workflow.

## 0. Expected input

The workflow uses sequencing data from the following SRA BioProjects:

* [PRJNA513137](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA513137)
* [PRJNA1072695](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1072695)

The workflow also requires a feature annotation file containing primer-template combinations, such as:

[`files/list_templates.txt`](files/list_templates.txt)


Expected columns:

```text
feature  template    primer
```


### 1. Preprocessing

File:

- [RMarkdown file](workflow/01_preprocessing_bash.Rmd)
- [HTML](https://clementeluangaspar.github.io/dePCR_drdorman/01_preprocessing_bash.html)

This file documents the preprocessing steps. I downloaded the SRA files, retrieved SRA/SRR metadata, converted SRA files to FASTQ, merged paired-end reads, and converted merged FASTQ files to FASTA/FNA format.

Main steps:

1. Download sequencing data from SRA.
2. Retrieve run metadata.
3. Convert SRA files to FASTQ.
4. Merge paired-end reads using VSEARCH.
5. Convert merged FASTQ files to FASTA/FNA.

### 2. Metrics, counts, and quality control

File:

```text
workflow/02_metrics_qc_counts.Rmd
```

This file summarizes preprocessing metrics and quality control results. It calculates the number of merged and unmerged reads, run FastQC and MultiQC, and ......

Main steps:

1. Run FastQC and MultiQC.
2. Summarize merged and unmerged reads.
3. Count primer-template features.
4. Generate the sample-by-feature count matrix.
5. Check basic count distributions and sample-level metrics.
6. .........

### 3. Data thinning

File:

```text
workflow/03_data_thinning.Rmd
```

This file applies data thinning to the count matrix. It splits the observed counts into conditionally related components that can support model fitting, validation, and/or selective inference strategies.

Main steps:

1. Load the primer-template count matrix.
2. ....

### 4. Negative binomial modeling with covariate-dependent overdispersion

File:

```text
workflow/04_nb_covariate_dependent_overdispersion.Rmd
```

This file fits negative binomial models in which the mean and the overdispersion parameter can depend on experimental covariates. I use these models to evaluate primer-template effects, experimental conditions, and overdispersion patterns.

Main steps:

1. Load the thinned count data.
2. .....
3. 
### 5. Figures and tables

File:

```text
workflow/05_figures_tables.Rmd
```

This file generates the final figures and tables for reporting. 

Main steps:

1. Load processed results.
2. ....

## License

Code in this repository is licensed under the MIT License.
