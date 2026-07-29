# dePCR Sequencing Workflow

This repository contains a reproducible workflow for processing dePCR sequencing data, generating primer-template count matrices, validating preprocessing results, and fitting Bayesian negative binomial models with covariate-dependent overdispersion.

## Repository description

The workflow starts with public sequencing data from NCBI SRA and proceeds through paired-end read merging, quality control, primer-template feature counting, metadata integration, data partitioning, Bayesian model fitting, model comparison, and out-of-sample validation.

The principal objective is to quantify how primer-template sequence configurations and mismatch characteristics affect observed read counts and count variability.

## Data sources

The workflow uses sequencing data from the following NCBI BioProjects:

- [PRJNA513137](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA513137)
- [PRJNA1072695](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1072695)

The analysis also uses:

- SRA run-information tables for both BioProjects;
- the supplementary sample metadata distributed with the original study;
- a primer-template feature annotation table;
- an experiment-specific feature table defining the expected features for each experiment;
- a feature metadata table containing mismatch status and nucleotide identity at the evaluated primer positions.

An example feature annotation file is provided at:

[`files/list_templates.txt`](files/list_templates.txt)

Expected columns:

```text
feature    template    primer
```

## Workflow organization

The analysis is documented in two main RMarkdown files.

### 1. Sequencing preprocessing

Files:

- [RMarkdown](workflow/01_preprocessing_bash.Rmd)
- [HTML](https://clementeluangaspar.github.io/dePCR_drdorman/01_preprocessing_bash.html)

This document describes the processing of raw sequencing data and generation of the sample-by-feature count matrix.

Main steps:

1. Download sequencing files from NCBI SRA.
2. Retrieve SRA and SRR run metadata.
3. Convert SRA files to paired FASTQ files.
4. Merge paired-end reads with VSEARCH.
5. Run FastQC and summarize the reports with MultiQC.
6. count merged and unmerged reads.
7. convert merged FASTQ files to FASTA/FNA.
8. count primer-template features in each sample.
9. generate the sample-by-feature count matrix.

The primer sequences are not trimmed before feature counting because primer identity is part of the feature definition. Each feature is defined by the simultaneous occurrence of the expected primer and template pattern in a merged read.

### 2. Count analysis and negative binomial modeling

Files:

- [RMarkdown](workflow/02_depcr_analysis.Rmd)
- [HTML](https://clementeluangaspar.github.io/dePCR_drdorman/02_depcr_analysis.html)

This document describes the statistical analysis performed after construction of the count matrix.

Main steps:

1. combine NCBI run metadata, locally observed read counts, and supplementary sample metadata;
2. verify that downloaded read totals agree with the published metadata;
3. calculate the total number and proportion of reads assigned to primer-template features;
4. reshape the count matrix from sample-wide to feature-long format;
5. retain the primer-template combinations expected for each experiment;
6. attach primer mismatch and nucleotide annotations;
7. construct sequence-group and mismatch-status covariates;
8. select experiment A10 and separate PCR and dePCR observations;
9. split biological or technical replicates into training and test sets;
10. fit Bayesian negative binomial models for the expected count;
11. allow the negative binomial shape parameter to depend on experimental covariates;
12. assess MCMC convergence and sampling quality;
13. compare supported models with approximate leave-one-out cross-validation;
14. evaluate the selected model on held-out replicates;
15. generate observed-versus-predicted figures and performance summaries.

## Statistical model

Counts are modeled using a negative binomial distribution with a log link.

For observation \(i\),

\[
Y_i \sim \mathrm{NegBinomial}(\mu_i, \phi_i),
\]

where \(\mu_i\) is the expected count and \(\phi_i\) is the negative binomial shape parameter. Under the parameterization used by `brms`,

\[
\mathrm{Var}(Y_i) = \mu_i + \frac{\mu_i^2}{\phi_i}.
\]

The mean model includes an offset for the total number of counted reads:

\[
\log(\mu_i) =
\log(\mathrm{total\ counted\ reads}_i)
+ \mathbf{x}_i^\top\boldsymbol{\beta}.
\]

The offset converts the model from an analysis of absolute counts into an analysis of relative feature abundance while retaining the original count-scale likelihood.

The shape parameter is modeled with a separate log-linear predictor:

\[
\log(\phi_i) = \mathbf{z}_i^\top\boldsymbol{\gamma}.
\]

A smaller value of \(\phi_i\) corresponds to greater extra-Poisson variation, whereas a larger value corresponds to lower overdispersion.

## Evaluated covariates

The current workflow evaluates the following predictors:

- `Temp`: amplification temperature, with 45°C as the reference level;
- `p3_status`, `mid_status`, and `p5_status`: match or mismatch status at the evaluated primer positions;
- `p3_base`, `mid_base`, and `p5_base`: nucleotide identity at each evaluated position;
- `group`: the combined nucleotide configuration `p3_base_mid_base_p5_base`.

The reference sequence group is:

```text
A_T_C
```

## Training and validation design

The current analysis focuses on experiment `A10`.

PCR and dePCR observations are separated before modeling. Six of the eight replicates are assigned to the training set, and the remaining two replicates are retained as an independent test set.

The same sampled replicate identifiers should be used for PCR and dePCR whenever the two methods are compared directly. This avoids creating method-specific train/test partitions.

## Candidate overdispersion models

The mean model uses the sequence group and, when appropriate, temperature:

```text
Counts ~ offset(log_total_counted) + Temp + group
```

Candidate shape models include:

```text
shape ~ 1
shape ~ Temp
shape ~ group
shape ~ p3_status
shape ~ mid_status
shape ~ p5_status
shape ~ p3_base
shape ~ mid_base
shape ~ p5_base
shape ~ p3_status + mid_status + p5_status
shape ~ Temp + p3_status
shape ~ Temp + mid_status
shape ~ Temp + p5_status
shape ~ Temp + p3_base
shape ~ Temp + mid_base
shape ~ Temp + p5_base
shape ~ Temp + group
shape ~ Temp + p3_status + mid_status + p5_status
```

Models are fitted with `brms` using the No-U-Turn Sampler. Model complexity is retained only when the chains converge and the additional terms improve predictive performance.

## Model diagnostics and comparison

A model is not interpreted solely because `brm()` completed successfully. The workflow examines:

- split-chain \(\widehat{R}\);
- bulk and tail effective sample sizes;
- divergent transitions;
- maximum tree-depth warnings;
- posterior predictive behavior;
- Pareto-\(k\) diagnostics from leave-one-out cross-validation.

Models with substantial convergence problems are excluded from biological interpretation and predictive comparison until they are reparameterized or refitted.

Among the supplied dePCR model runs, the simpler models with constant shape or shape depending individually on mismatch status showed substantially better convergence than several highly parameterized shape models. Models with `shape ~ group`, `shape ~ Temp + group`, and several PCR shape models displayed high \(\widehat{R}\) and very small effective sample sizes and therefore should not be interpreted in their current form.

## Held-out test evaluation

The selected training model is evaluated on the two held-out replicates using:

- posterior expected counts;
- posterior predictive intervals;
- root mean squared error;
- mean absolute error;
- Pearson correlation;
- Spearman correlation;
- empirical coverage of the 95% posterior predictive interval;
- total and mean log predictive density.

Observed-versus-predicted plots are generated on both the original and logarithmic scales.

## Repository structure

```text
dePCR_drdorman/
├── README.md
├── LICENSE
├── files/
│   ├── list_templates.txt
│   ├── dePCR_feature_metadata.xlsx
│   ├── dePCR_experiment_features.txt
│   └── other required metadata files
├── workflow/
│   ├── 01_preprocessing_bash.Rmd
│   └── 02_depcr_analysis.Rmd
├── results/
│   ├── metadata/
│   ├── model_fits/
│   ├── model_summaries/
│   ├── predictions/
│   ├── tables/
│   └── figures/
└── docs/
    ├── 01_preprocessing_bash.html
    └── 02_depcr_analysis.html
```

Large raw FASTQ, FNA, and fitted-model files do not need to be committed to GitHub. The repository should instead preserve the scripts, small metadata files, summaries, figures, and instructions needed to recreate them.

## Main preprocessing outputs

The preprocessing workflow generates:

```text
sra_runinfo_PRJNA513137.csv
sra_runinfo_PRJNA1072695.csv
fastq/
merged_vsearch/
merge_logs/
fastqc_reports/
merge_summary.txt
MergedFASTAFiles/
FinalCount_samples_rows.tsv
```

## Main analysis outputs

The analysis workflow generates:

```text
metadata_read_comparison.tsv
sample_feature_counts_long.tsv
train_replicate_ids.txt
model_diagnostics.tsv
loo_model_comparison.tsv
test_predictions.tsv
test_performance_metrics.tsv
observed_vs_predicted_test.png
observed_vs_predicted_by_group.png
```

Bayesian fit objects and text summaries may be stored separately:

```text
results/model_fits/*.rds
results/model_summaries/*_summary.txt
```

## Software

The workflow uses:

### Command-line tools

- NCBI SRA Toolkit
- VSEARCH
- FastQC
- MultiQC
- seqtk
- standard Bash utilities

### R packages

- `readxl`
- `ShortRead`
- `dplyr`
- `tidyr`
- `ggplot2`
- `brms`
- `posterior`
- `loo`
- `matrixStats`

## Reproducibility notes

File paths in the analysis document are defined in a single setup section and should be adapted to the local project directory.

The random seed used for replicate partitioning and Bayesian fitting is:

```r
10231991
```

The workflow saves the selected training-replicate identifiers so that the same split can be reused in later analyses.

Computationally intensive model-fitting chunks are controlled with the RMarkdown parameter:

```yaml
params:
  run_models: false
```

Set `run_models: true` when the Bayesian models need to be fitted again. When it is `false`, the report reads previously saved model objects and summaries.

## License

Code in this repository is licensed under the MIT License.
