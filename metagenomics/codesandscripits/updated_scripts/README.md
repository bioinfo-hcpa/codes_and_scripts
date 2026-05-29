# Metagenomic Analysis Workflows

This directory contains standardized workflows for the analysis of 16S rRNA amplicon sequencing data developed by the HCPA Bioinformatics Center.

The scripts are organized according to the main stages of microbiome data analysis, from sequence preprocessing to downstream statistical analyses for both cross-sectional and longitudinal study designs.

## Directory Structure

### [preprocessing_scripts](preprocessing_scripts)

Workflows for preprocessing raw sequencing data using **DADA2** and **phyloseq**.

Main objectives:

* Quality assessment and filtering of sequencing reads.
* Sequence denoising and ASV inference.
* Taxonomic assignment.
* Generation of analysis-ready phyloseq objects.

Available workflows:

* Single sequencing run datasets.
* Multi-run sequencing datasets (BIGDATA).
* Ion Torrent sequencing datasets.

---

### [crossSectional_scripts](crossSectional_scripts)

Workflows for downstream microbiome analyses in cross-sectional studies.

Main analyses:

* Taxonomic composition overview.
* Relative abundance profiling.
* Alpha and Beta diversity analysis.
* Differential abundance testing using ANCOM-BC2.

Main statistical approaches:

* Shannon diversity index.
* Mann-Whitney/ Kruskal–Wallis tests.
* PCoA and NMDS ordinations.
* β-dispersion analysis.
* PERMANOVA.
* Differential abundance analysis using ANCOM-BC2.

---

### [timeSeries_scripts](timeSeries_scripts)

Workflows for downstream microbiome analyses in longitudinal (time-series) studies.

Main analyses:

* Temporal taxonomic characterization.
* Alpha and Beta diversity analysis.
* Differential abundance testing using MaAsLin2.

Main statistical approaches:

* Shannon diversity index.
* Linear Mixed-Effects Models (LME).
* Continuous autoregressive correlation structures (CAR(1)).
* PCoA and NMDS ordinations.
* Longitudinal modeling of ordination axes.
* β-dispersion analysis.
* PERMANOVA.
* Temporal autocorrelation analyses.
* Multivariable association analysis.
* Differential abundance analysis using MaAsLin2.

---

## RECOMMEND WORKFLOW

1. Run the appropriate preprocessing workflow in `preprocessing_scripts`.
2. Generate the processed phyloseq object.
3. Select the downstream workflow according to the study design:
   * `crossSectional_scripts` for cross-sectional studies.
   * `timeSeries_scripts` for longitudinal/ time series studies.
4. Perform diversity analyses and differential abundance testing.
5. Generate publication-ready tables and figures.

## Main R Packages

* `DADA2`
* `phyloseq`
* `microbiome`
* `vegan`
* `ANCOM-BC2`
* `MaAsLin2`
* `compositions`
* `dplyr`
* `tidyverse`
* `ggplot2`

---

## Maintenance and Updates

- Original workflows: **HCPA Bioinformatics Center**
- Reviewer: **Thaiane Nascimento**, Bioinformatician (Key workflow updates and documentation)
- Last update: **March 2026**

