# Metagenomic Cross-Sectional Analysis Scripts

This directory contains workflows for downstream analysis of 16S rRNA amplicon sequencing data from cross-sectional studies, developed by the HCPA Bioinformatics Center.
These workflows are intended for datasets previously processed using the preprocessing pipelines available in the *preprocessing_scripts* directory.

## Analysis Overview

The analyses are primarily performed using the **phyloseq**, **microbiome**, **vegan**, **ANCOMBC**, **dplyr**, **tidyverse**, and **ggplot2** R packages. The workflows provide standardized approaches for microbiome characterization, diversity analyses, and differential abundance testing.

## Available Scripts

### `2_HCPA_NBioinfo_16S_crossSectional_overview.Rmd`

Descriptive analysis workflow for microbiome characterization across experimental groups.

Main objectives:

* Characterize the overall microbial composition of the dataset.
* Remove low-prevalence and low-abundance taxa to prioritize biologically relevant microorganisms.
* Generate relative abundance summaries at the phylum, family, and genus levels.
* Produce descriptive tables and visualizations for group comparisons.

Main packages:
* `phyloseq`
* `dplyr`
* `tidyverse`
* `ggplot2`

Outputs:
* Relative abundance tables.
* Stacked bar plots.
* Taxonomic composition summaries.

---
### `3_HCPA_NBioinfo_16S_crossSectional_process_alphaBetadiv.Rmd`

Workflow for alpha and beta diversity analyses.

#### Alpha Diversity

Alpha diversity evaluates microbial diversity within individual samples.

Methods:
* Shannon diversity index.
* Group comparisons using the non-parametric Kruskal–Wallis test.

Outputs:
* Diversity tables.
* Boxplots with statistical significance annotations.

#### Beta Diversity

Beta diversity evaluates differences in microbial community composition between samples and experimental groups.

Methods:
* Compositional approaches using the `compositions` package.
* Ecological distance-based analyses using the `vegan` package.
* Ordination techniques:
  * PCA
  * PCoA
  * NMDS

Statistical testing:
* β-dispersion analysis (`betadisper`) followed by ANOVA.
* PERMANOVA (`adonis2`) for global community comparisons.

Outputs:
* PCA, PCoA, and NMDS plots.
* 95% confidence ellipses.
* Explained variance (PCA and PCoA).
* Stress values (NMDS).
* Exported sample coordinates.
* Statistical summary tables containing F statistics, R² values, and p-values.

---

### `4_HCPA_NBioinfo_16S_crossSectional_diffAbundance.Rmd`

Workflow for differential abundance analysis.

Main objective:
* Identify microbial genera with significantly different abundances between experimental groups.

Methods:
* ANCOM-BC2 implemented through the `ANCOMBC` package.
* Genus-level analysis.
* Correction for compositional bias and library size differences.
* Filtering of rare taxa.
* Multiple testing correction.
* Statistical significance threshold of 5%.

Outputs:
* Tables of differentially abundant genera.
* Direction of abundance change between groups.
* Adjusted p-values.
* Horizontal bar plots summarizing differential abundance results.
