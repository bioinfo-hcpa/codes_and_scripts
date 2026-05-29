# Metagenomic Time-Series Analysis Scripts

This directory contains workflows for downstream analysis of 16S rRNA amplicon sequencing data from longitudinal (time-series) studies, developed by the HCPA Bioinformatics Center.

These workflows are intended for datasets previously processed using the *Metagenomic Preprocessing Scripts* directory.

## Analysis Overview

The analyses are primarily performed using the **phyloseq**, **microbiome**, **vegan**, **MaAsLin2**, **dplyr**, **tidyverse**, and **ggplot2** R packages. The workflows provide standardized approaches for microbiome characterization, diversity analyses, and differential abundance testing across multiple time points.

## Available Scripts

### [2_HCPA_NBioinfo_16S_timeSeries_overview.Rmd](2_HCPA_NBioinfo_16S_timeSeries_overview.Rmd)

Descriptive analysis workflow for microbiome characterization across study time points and experimental groups.

Main objectives:
* Characterize the overall microbial composition throughout the study period.
* Remove low-prevalence and low-abundance taxa to prioritize biologically relevant microorganisms.
* Generate relative abundance summaries at the phylum, family, and genus levels.
* Explore temporal patterns in microbial composition.
* Produce descriptive tables and visualizations for comparisons across time points and study groups.

Main packages:
* `phyloseq`
* `dplyr`
* `tidyverse`
* `ggplot2`

Outputs:
* Relative abundance tables.
* Stacked bar plots.
* Taxonomic composition summaries.
* Temporal abundance visualizations.

---

### [3_HCPA_NBioinfo_16S_timeSeries_process_alphaBetadiv.Rmd](3_HCPA_NBioinfo_16S_timeSeries_process_alphaBetadiv.Rmd)

Workflow for alpha and beta diversity analyses in longitudinal microbiome studies.

#### Alpha Diversity

Alpha diversity evaluates microbial diversity within individual samples over time.

Methods:
* Shannon diversity index.
* Longitudinal comparisons across time points and experimental groups.
* Linear Mixed-Effects Models (LME) including:
  * Fixed effects for Timepoint, Group, and Timepoint vs Group interaction.
  * Random intercepts for individual animals.
  * Continuous autoregressive correlation structure (CAR(1)) to account for temporal autocorrelation among repeated measurements.

Additional analyses:
* Residual diagnostics using the Shapiro-Wilk test.
* Confidence intervals for model coefficients.
* Predicted marginal effects visualization.
* Quantification of between-animal variability.
* Evaluation of temporal autocorrelation strength.

Outputs:
* Shannon diversity tables.
* Mixed-effects model summaries.
* Confidence intervals of model coefficients.
* Predicted marginal effects plots.
* Residual diagnostics reports.
* Temporal autocorrelation estimates.
* Publication-ready boxplots and longitudinal visualizations.

#### Beta Diversity

Beta diversity evaluates longitudinal changes in microbial community composition across samples, groups, and time points.

Methods:
* Aitchison distance for compositional analyses.
* Jaccard distance for presence/absence analyses.
* Ordination techniques:
  * PCoA
  * NMDS

Community-level statistical testing:
* β-dispersion analysis (`betadisper`) followed by ANOVA.
* PERMANOVA (`adonis2`) for global microbiome comparisons.

Longitudinal modeling:

* Linear Mixed-Effects Models applied to ordination axes.
* Fixed effects: Timepoint, Group, Timepoint vs Group interaction;
* Random intercepts for individual animals.
* Continuous autoregressive correlation structure (CAR(1)) to model temporal dependence.

Additional analyses:
* Axis-specific interpretation of temporal microbiome trajectories.
* Quantification of between-animal variability.
* Evaluation of within-animal temporal autocorrelation.
* Spaghetti plots for ordination dynamics over time.

Outputs:
* PCoA and NMDS plots.
* 95% confidence ellipses.
* Explained variance (PCoA).
* Stress values (NMDS).
* Exported ordination coordinates.
* Mixed-effects model summaries for ordination axes.
* Temporal autocorrelation statistics.
* Spaghetti plots of longitudinal microbiome trajectories.
* Statistical summary tables containing estimates, standard errors, F statistics, R² values, and p-values.


---

### [4_HCPA_NBioinfo_16S_timeSeries_diffAbundance.Rmd](4_HCPA_NBioinfo_16S_timeSeries_diffAbundance.Rmd)

Workflow for differential abundance analysis in longitudinal microbiome studies.
* Identify microbial genera significantly associated with study time points, experimental groups, or other metadata variables.

Methods:
* Differential abundance analysis using the `MaAsLin2` package.
* Multivariable association modeling.
* Genus-level analysis.
* Inclusion of fixed and/or random effects when appropriate.
* Filtering of low-abundance taxa.
* Multiple testing correction using false discovery rate (FDR).
* Statistical significance threshold of 5%.

Outputs:
* Tables containing significantly associated microbial genera.
* Effect sizes, coefficients, and adjusted p-values (q-values).
* Summary tables of model results.
* Publication-ready visualizations of significant associations.

---

## Maintenance and Updates

- Original workflows: **HCPA Bioinformatics Center**
- Reviewer: **Thaiane Nascimento**, Bioinformatician (Key workflow updates and documentation)
- Last update: **March 2026**

