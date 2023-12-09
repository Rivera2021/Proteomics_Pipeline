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
    - Params:
        - Metadata_path: Path to Metadata in .xlsx format
        - Count_path: Path to matrix with reads count in .tsv format
        - Match_feature: Column name from metadata that should coincide with the sample names used in the count matrix
        - Output_file_path: Path to a directory to store results
        - SampleNamepath: Path to .csv file containing the metadata column names that will be used to uniquely identified a sample using biological features
    - Output:
        - Import_Data.xlsx: File with Count matrix and metadata after a few modifications needed to start the pipeline
* QC_PRENORMALIZATION: QC based on the quality of the sequencing
    - Params:
          - Output_file_path: Path to a directory to store results
          - nMust: Number of IQRs below Q1(25% quantile) to consider outlier in number of reads and number of non-zero genes. 2 recommended 
    - Output:
        - Outliers.xlsx: List of samples selected as outliers based on the QC sequencing depth and number of non-zero genes
        - QC_Data.xlsx: File with Count Matrix and Metadata without outliers selected
        - QC_Prenormalization: QC plots
* PRE_FILTERING: Filtering genes depending on different criteria 
    - Params:
        - Output_file_path: Path to a directory to store results
        - PRE_FILTER: PREV (Prevalence of gene across samples is used as a prefiltering criterion)
        - Prev_perc: Fraction of samples in which a gene has to be detected to be considered as part of the analysis
    - Output:
        - Prefilter_Data.xlsx: File with Count matrix with gene prefiltering applied and Metadata.
* DESEQ_NORM: DESeq normalization and PCA plots
    - Params:
        - Output_file_path: Path to a directory to store results
        - QCNORM: Whether the normalization is made before or after the QC based on gene expression behavior. PRE_QCNORM (Normalization before QC), POST_QCNORM (Normalization after QC)
        - CoarseConditions: Features from metadata to build design matrix of DESeq normalization. These features will be used to select comparisons
        - METHOD_NORM: Normalization method. Standard (Normalization obtained from DESeq object)
        - VST_FILTER: Whether to apply Variance Stabilizing transformation to the normalized matrix. VST_ON or VST_OFF. VST_ON recommended
        - SVD_FILTER: Whether to apply SVD truncation to eliminate noisy PCA direction contributions. SVD_OFF or SVD_ON
        - PlotPCA: Whether to produce PCA plots. 'PCA_PLOT'
        - Control_Neg_PCA: Negative control to show in PCAs. Use "" if only samples associated to a particular treatment are desired in the plot
        - Control_Pos_PCA: Positive control to show in PCA. Use "" if there is none.
    - Output:
      
        - DESeq_Norm.RData: DESeq2 object with normalization and DEG information
        - Norm_Data.xlsx: Count matrix after applying vst on the DESeq2 normalized counts
        - PCA_PLOTS: PCA plots using Norm_Data per treatment and for all samples
          
* QC_POSTNORMALIZATION: Outlier detection based on gene expression and within replicate correlation
    - Params:
      
        - Output_file_path: Path to a directory to store results
        - OUTLIER_FILTER: Method to detect outliers. GENTLE_REP (Method based on correlation with replicates only)
        - RepFeatures_path: Metadata features to identify biological replicates 
   
    - Output:
      
        - Outlier_DF.xlsx:
            - Bad_Samples_cond: Samples associated to conditions where correlation among replicates is significantly low
            - Outlier_samples: Samples that correlate very poorly with all its replicates
            - Outlier_Selected: Bad_Samples + Outlier_samples. Samples to be removed from the rest of the analysis.
        - Correlation_within_replicates.pdf: Plots showing correlation within replicates across all conditions
        - PCA_WITH_OUTLIERS: PCA plots showing selected outliers
  
* DEG_FUNCTION: Finds differentially expressed genes

    - Output: A folder per molecule
        - DEG_MXXX_.RData: List containing DEG for every comparison containing the selected molecule
        - Volcano plots for every comparison
          
  
    

