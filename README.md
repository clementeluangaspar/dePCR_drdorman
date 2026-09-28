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

~~~text
feature    template    primer
~~~

## Workflow organization

The analysis is documented in two main RMarkdown files and a simplified R script describing the main statistical workflow.

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
5. Run FastQC and summarize quality-control results.
6. Count merged and unmerged reads.
7. Convert merged FASTQ files to FASTA/FNA.
8. Count primer-template features in each sample.
9. Generate the sample-by-feature count matrix.

The primer sequences are not trimmed before feature counting because primer identity is part of the feature definition. Each feature is defined by the simultaneous occurrence of the expected primer and template pattern in a merged read.

### 2. Statistical analysis

Files:

- [RMarkdown](workflow/02_depcr_analysis.Rmd)
- [HTML](https://clementeluangaspar.github.io/dePCR_drdorman/02_depcr_analysis.html)
- [`02_modeling.R`](workflow/02_modeling.R)

This part of the workflow describes the exploratory analysis and Bayesian modeling performed after construction of the count matrix.

Main steps:

1. Combine sequencing counts with experiment and feature metadata.
2. Calculate assigned-read proportions and exploratory summaries.
3. Restrict primary model development to experiments A10 and A64.
4. Analyze PCR and DePCR separately.
5. Remove sequencing libraries with 13,000 or fewer assigned reads.
6. Construct alternative representations of primer sequence information.
7. Fit Bayesian distributional negative binomial candidate models.
8. Compare models using PSIS-LOO.
9. Fit the selected final model for each protocol.
10. Compare the final distributional negative binomial model with a conventional negative binomial model using a constant shape parameter.
11. Apply the fitted models to A27, B10, and B27 for exploratory transfer assessment.
12. Generate posterior summaries and observed-versus-predicted figures.

## Experimental datasets

The primary modeling analysis uses experiments A10 and A64:

- `A10`: ST0 template and a 10-primer pool;
- `A64`: ST0 template and all 64 primer variants.

Additional experiments are used for exploratory transfer assessment:

- `A27`: ST0 template and an alternative 27-primer pool;
- `B10`: ST0-ST9 template mixture and a 10-primer pool;
- `B27`: ST0-ST9 template mixture and a 27-primer pool.

PCR and DePCR are modeled separately because locus-specific primers play different roles in the two amplification protocols.

## Statistical preprocessing

The primary modeling dataset is restricted to A10 and A64.

Sequencing libraries with 13,000 or fewer reads assigned to the modeled primer-template features are excluded.

After filtering, the final modeling datasets contain:

- 1,184 observations for PCR;
- 898 observations for DePCR.

Sequencing depth is accounted for using an offset rather than rarefying libraries to a common depth.

## Statistical model

Feature counts are modeled using a negative binomial distribution:

$$
Y_i \sim \mathrm{NegBinomial}(\mu_i,\phi_i),
$$

where $\mu_i$ is the expected count and $\phi_i$ is the negative binomial shape parameter.

Under the parameterization used by `brms`,

$$
\mathrm{Var}(Y_i)
=
\mu_i+\frac{\mu_i^2}{\phi_i}.
$$

For a fixed mean, smaller values of $\phi_i$ correspond to greater extra-Poisson variation, whereas larger values indicate lower residual overdispersion.

The mean model includes an offset for sequencing depth:

$$
\log(\mu_i)
=
\log(L_i)
+
\mathbf{x}_i^\top\boldsymbol{\beta},
$$

where $L_i$ is the total number of reads assigned to the modeled primer-template features for the corresponding sequencing library.

The shape parameter is modeled separately:

$$
\log(\phi_i)
=
\mathbf{z}_i^\top\boldsymbol{\gamma}.
$$

This distributional specification allows predictors of expected abundance and predictors of residual variability to differ.

## Primer sequence representations

Three position-specific representations are evaluated:

- nucleotide identity at each controlled position;
- mutation type at each controlled position;
- match-versus-mismatch status at each controlled position.

For each representation, increasing amounts of positional information are considered:

- 3' position only;
- 3' and middle positions;
- 3', middle, and 5' positions.

Each representation is evaluated in:

- the mean model;
- the shape model;
- both the mean and shape models.

In addition, the complete three-position nucleotide configuration is represented as a single joint primer state.

Because each position can contain one of four nucleotides, the complete joint representation contains up to:

~~~text
4 x 4 x 4 = 64
~~~

primer configurations.

The reference configuration is:

~~~text
A_T_C
~~~

A total of 30 candidate models are evaluated separately for PCR and DePCR.

## Evaluated covariates

The workflow evaluates the following variables:

- `Exp`: experiment, with A10 as the reference level;
- `Temp`: annealing temperature, with 45°C as the reference level;
- `p3_status`: match or mismatch status at the 3' position;
- `mid_status`: match or mismatch status at the middle position;
- `p5_status`: match or mismatch status at the 5' position;
- `p3_base`: nucleotide identity at the 3' position;
- `mid_base`: nucleotide identity at the middle position;
- `p5_base`: nucleotide identity at the 5' position;
- mutation-specific variables for the corresponding positions;
- `group`: the complete three-position nucleotide configuration.

## Candidate-model comparison

Candidate models are designed to evaluate whether primer information improves predictive performance by explaining:

1. expected abundance;
2. residual variability;
3. both components simultaneously.

Experiment and annealing temperature are included as experimental covariates, and the mean model includes the sequencing-depth offset.

Candidate models are compared using Pareto-smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO).

Model comparisons are based on differences in expected log predictive density (ELPD).

More negative values of $\Delta$ELPD indicate worse predictive performance relative to the best-performing model.

## Final model

The same final model structure is used for PCR and DePCR, but the two protocols are fitted separately.

The final mean model is:

~~~text
Counts ~ offset(log_total_counted) + Exp + Temp + group
~~~

The final shape model is:

~~~text
shape ~ Exp + Temp + p3_status + mid_status + p5_status
~~~

Reference levels are:

- experiment: `A10`;
- annealing temperature: `45C`;
- mismatch status: `Match`;
- primer configuration: `A_T_C`.

## Bayesian fitting

Models are fitted in `brms` using the negative binomial family and the No-U-Turn Sampler through `rstan`.

Final models use:

~~~text
chains = 4
iter = 4000
warmup = 1000
post-warmup draws = 12000
seed = 10231991
adapt_delta = 0.99
max_treedepth = 15
~~~

The main software versions used for the final analysis are:

~~~text
brms        2.23.0
rstan       2.32.7
StanHeaders 2.39.1
~~~

Posterior summaries are reported using posterior means and 95% credible intervals.

## Model diagnostics

Sampling quality is evaluated using:

- split R-hat;
- bulk effective sample size;
- tail effective sample size;
- Pareto-k diagnostics from PSIS-LOO.

Models with substantial convergence problems are not used for biological interpretation or final predictive comparison.

## Constant-shape negative binomial comparison

For both PCR and DePCR, the final distributional negative binomial model is compared with a conventional negative binomial model using the same mean structure but a constant shape parameter.

The conventional negative binomial model uses:

~~~text
Counts ~ offset(log_total_counted) + Exp + Temp + group
shape ~ 1
~~~

The distributional model uses:

~~~text
Counts ~ offset(log_total_counted) + Exp + Temp + group
shape ~ Exp + Temp + p3_status + mid_status + p5_status
~~~

The models are compared using PSIS-LOO to evaluate whether allowing residual overdispersion to vary across experimental conditions improves predictive performance.

## Exploratory transfer assessment

The final models are developed using A10 and A64 and are subsequently applied to A27, B10, and B27.

A27 retains the ST0 template but uses a different primer pool.

B10 and B27 contain mixtures of templates ST0-ST9 and therefore introduce template-specific sequence context and competition not explicitly represented in the A10/A64 development system.

For visualization, observed counts and posterior expected counts are converted to within-template proportions.

For a given template,

$$
p_k
=
\frac{\mathrm{Count}_k}
{\sum_j \mathrm{Count}_j}.
$$

Posterior expected counts are normalized in the same way.

Observed and predicted proportions are compared using an identity line.

The transfer analysis is exploratory rather than a formal external validation and is used to identify template-, temperature-, and mismatch-dependent departures from the relationships estimated in the A10/A64 system.

## Posterior prediction

Posterior expected counts are obtained using:

~~~r
posterior_epred()
~~~

When posterior predictive observations are required, they are generated using:

~~~r
posterior_predict()
~~~

Posterior expected values are used for the main observed-versus-predicted transfer figures.

## Main preprocessing outputs

The preprocessing workflow generates files including:

~~~text
sra_runinfo_PRJNA513137.csv
sra_runinfo_PRJNA1072695.csv
fastq/
merged_vsearch/
merge_logs/
fastqc_reports/
merge_summary.txt
MergedFASTAFiles/
FinalCount_samples_rows.tsv
~~~

## Main analysis outputs

The statistical workflow produces objects and summaries including:

~~~text
model_comparison_PCR.tsv
model_comparison_DePCR.tsv
component_gain_PCR.tsv
component_gain_DePCR.tsv
final_model_PCR.rds
final_model_DePCR.rds
prediction_A27.tsv
prediction_B10.tsv
prediction_B27.tsv
~~~

Model objects and summaries can be stored under:

~~~text
results/model_fits/
results/model_summaries/
~~~

Figures generated by the analysis include:

~~~text
model_comparison_PCR.png
model_comparison_DePCR.png
component_gain_PCR.png
component_gain_DePCR.png
pcr_prediction_a27_temp_status.png
depcr_prediction_a27_temp_status.png
pcr_prediction_b10_temp_status.png
depcr_prediction_b10_temp_status.png
pcr_prediction_b27_temp_status.png
depcr_prediction_b27_temp_status.png
~~~

## Simplified modeling script

A simplified R script is included at:

[`workflow/02_modeling.R`](workflow/02_modeling.R)

The purpose of this script is to document the principal statistical workflow without reproducing every exploratory model and plotting step contained in the full RMarkdown analysis.

The script includes examples of:

- preparation of PCR and DePCR datasets;
- definition of reference levels;
- final `brms` model formulas;
- distributional negative binomial fitting;
- constant-shape negative binomial fitting;
- PSIS-LOO comparison;
- posterior summaries;
- posterior expected-value calculation.

The complete candidate-model analysis, exploratory summaries, and figure-generation workflow are documented in the RMarkdown files.

## Reproducibility

The repository is intended to provide a reproducible record of the major preprocessing and statistical steps used in the analysis.

Sequencing preprocessing is documented in:

~~~text
workflow/01_preprocessing_bash.Rmd
~~~

Statistical analysis is documented in:

~~~text
workflow/02_depcr_analysis.Rmd
~~~

A shorter reference implementation of the main Bayesian modeling steps is provided in:

~~~text
workflow/02_modeling.R
~~~

Because Bayesian model fitting can be computationally intensive, fitted model objects may be stored separately from the rendered documentation.

## License

Code in this repository is licensed under the MIT License.
