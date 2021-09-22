[![DOI](https://zenodo.org/badge/366495264.svg)](https://zenodo.org/badge/latestdoi/366495264)

# COVID_ML

This repository contains the de-identified data and code used in our retrospective study: "Unbiased identification of clinical characteristics predictive of COVID-19 severity"

https://pubmed.ncbi.nlm.nih.gov/34089403/

## Installation guide
No installation necessary

## Demo
Example de-identified data is included in the `deID.clem.xlsx` file.  Note that some variables have been removed to preserve patient privacy.

## Instructions
1. Download `deID.clem.xlsx` and `deID.analysis.R` files
2. Run `deID.analysis.R` in Rstudio.  This script will generate the figures presented in our manuscript.  Note that script has been modified to account for variables removed from dataset for patient privacy.

## System requirements and software
R	The Comprehensive R Archive Network	v3.6.1\
Rstudio 	Rstudio, Inc.	v1.2.1578\
ggpubr (version 0.4.0)\
corrplot (version 0.84)\
survminer (version 0.4.8) \
survival (version 3.2-7)\
ggplot2 (version 3.3.0)\
pheatmap (version 1.0.12)\
EnhancedVolcano (version 1.4.0) 
