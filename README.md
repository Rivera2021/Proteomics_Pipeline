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
* raw: Contains the raw count matrix and metadata to start the analysis


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
        - PCA_PLOTS: PCA plots using Norm_Data per treatment and for all samples
* QC_POSTNORMALIZATION: Outlier detection based on gene expression and within replicate correlation

    - Output:
        - Outlier_DF.xlsx: MH_dist_df (Mahalanobis distance with respect to the cloud for every sample), MH_crt (Mahalanobis distance threshold corresponding to P < 0.001), Corr (Correlation matrix across all samples), Outlier_Over (Outliers list first pass), Outlier_RepV2 (Outliers list second pass, based only on correlation between replicates), Outlier_Selected (All outliers selected to be removed)
        - Outliers_Selected.csv: All outliers selected to be removed
        - Mahalanobis_screeplot: Rnak plot of MH_dist_df
        - MAHALANOBIS_OUTLIERS_PCA: PCA showing outliers identified using Mahalanobis distance. Note: Not all these samples will be removed in the GENTLE model
        - SELECTED_OUTLIERS_PCA: PCA showing outliers removed from the rest of the analysis
  
* DEG_FUNCTION: Finds differentially expressed genes

    - Output: A folder per molecule
        - DEG_MXXX_.RData: List containing DEG for every comparison containing the selected molecule
        - Volcano plots for every comparison
          
  
    

