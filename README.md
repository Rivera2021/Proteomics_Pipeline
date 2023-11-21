# Transcriptomics_Pipeline

This project is used to create the Nextflow pipeline for bulk RNAseq transcriptomics data analysis. 

Started on 2023-10-27.

## Directory structure

* **R** - Resuable R code (functions etc.)
* **analysis** - R Markdown analysis files
* **docs** - Rendered analysis reports
* **data** - Raw data used for analysis

## Data

* Data_for_pipeline: Contains different type of data that is used along the pipeline
* raw: Contains the raw count matrix and metadata to start analysis


## Analysis


## Code

* IMPORT_DATA.R: Validates that Count matrix and Metadata agree with each other
* QC_PRENORMALIZATION: QC based on the quality of the sequencing
* PRE_FILTERING: Filtering genes depending on different criteria
* DESEQ_NORM: DESeq normalization and PCA plots
* QC_POSTNORMALIZATION: Outlier detection based on gene expression and within replicate correlation
* DEG_FUNCTION: Finds differentially expressed genes

