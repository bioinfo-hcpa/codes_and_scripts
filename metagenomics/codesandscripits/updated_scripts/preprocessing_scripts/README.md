# Metagenomic Preprocessing Scripts

This directory contains preprocessing workflows for 16S rRNA amplicon sequencing data developed by the HCPA Bioinformatics Center.

The workflows are primarily based on the DADA2 and phyloseq R packages and provide standardized preprocessing steps prior to downstream microbiome analyses, including quality assessment, sequence denoising, taxonomic assignment, and generation of analysis-ready datasets.

## Available Scripts

- ### [1-HCPA_NBioinfo_16S_preprocess_singleRun.Rmd](1-HCPA_NBioinfo_16S_preprocess_singleRun.Rmd)
Workflow for datasets generated from a single sequencing run, where all sequencing files are located within a single directory.

- ### [1-HCPA_NBioinfo_16S_preprocess_BIGDATA.Rmd](1-HCPA_NBioinfo_16S_preprocess_BIGDATA.Rmd)        
Workflow for datasets generated across multiple sequencing runs, where sequencing files are distributed across two or more directories.
The script is designed to identify, import, and process samples from multiple locations within a unified workflow.

- ### [1-HCPA_NBioinfo_16S_preprocess_ionTorrent.Rmd](1-HCPA_NBioinfo_16S_preprocess_ionTorrent.Rmd)
Workflow for preprocessing 16S rRNA sequencing data generated using the Ion Torrent platform.

## Maintenance and Updates

- Original workflows: **HCPA Bioinformatics Center**
- Reviewer: **Thaiane Nascimento**, Bioinformatician (Key workflow updates and documentation)
- Last update: **March 2026**
