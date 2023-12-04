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

* IMPORT_DATA.R: Validates that Count matrix and Metadata agree with each other.
  
    - Output:
        - Import_Data.xlsx: File with Count matrix and metadata after a few modifications needed to start the pipeline
* QC_PRENORMALIZATION: QC based on the quality of the sequencing

    - Output:
        - Outliers.xlsx: List of samples selected as outliers based on the QC sequencing depth and number of non-zero genes
        - QC_Data.xlsx: File with Count Matrix and Metadata without outliers selected
        - QC_Prenormalization: QC plots
* PRE_FILTERING: Filtering genes depending on different criteria 

    - Output:
        - Prefilter_Data.xlsx: File with Count matrix with gene prefiltering applied and Metadata.
* DESEQ_NORM: DESeq normalization and PCA plots

    - Output:
        - DESeq_Norm.RData: DESeq2 object with normalization and DEG information
        - Norm_Data.xlsx: Count matrix after applying vst on the DESeq2 normalized counts
        - PCA_PLOTS: PCA plots using Norm_Data for all samples and per treatment
* QC_POSTNORMALIZATION: Outlier detection based on gene expression and within replicate correlation

    - Output:
        - DEG_MOLXXX__.RData: 
  
* DEG_FUNCTION: Finds differentially expressed genes

