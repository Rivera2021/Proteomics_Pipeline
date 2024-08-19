# Transcriptomics_Pipeline

This project is used to create the Nextflow pipeline for bulk RNAseq transcriptomics data analysis. 

Started on 2023-10-27.

## Directory structure

* **R** - Resuable R code (functions etc.)
* **data** - Raw data used for analysis

## Data

* Data_for_pipeline: Contains different type of data that is used along the pipeline
* raw: Contains the raw count matrix and metadata to start the analysis


## Analysis


## Code

Three reports are generated:
* QC_REPORT: Contains technical QC, PCA plots, and correlations plots. Uses the following files
    - QC_REPORT.Rmd
        -Expects to find Count matrix and Metadata in an xlsx format within IMPORT_DATA folder (This will should be optimize once protocol input formats and metadata have been totally determined). Count matrix requires a column named 'Gene'
        
    - QC_PRENORMALIZATION_CL
    - PRE_FILTERING
    - DESEQ_NORM
 
* DEG_ENRICH: For every cell line and every molecule a report is generated. It contains DEGs and GSEA enriched pathways. Uses the following files
    - MAIN_DEG_ENRICH_perdrug.R
    - DEG_ENRICH_perdrug.Rmd
    - DEG_FUNCTION_DA.R
    
* MAIN_INTEGRATION: Contains integration with chemoproteomics across selected conditions. Uses the following files
    - MAIN_INTEGRATION.R
    - INTEGRATION_REPORT.Rmd

The following functions have been used throughout the pipeline


* IMPORT_DATA.R: Validates consistency between count matrix and Metadata
    - Params:
        - Metadata_path: Path to Metadata in .xlsx format
        - Count_path: Path to matrix with read count in .tsv format
        - Match_feature: Column name from metadata that matched the sample names used in the count matrix
        - Output_file_path: Directory path for storing results
        - SampleNamepath: Path to .csv file containing the metadata column names that will be used to uniquely identify a sample using biological features
    - Output:
        - Import_Data.xlsx: File with count matrix and metadata after a few modifications needed to start the pipeline
* QC_PRENORMALIZATION: QC based on the quality of the sequencing
    - Params:
        - Output_file_path: Directory path for storing results
        - nMust: Number of IQRs below Q1(25% quantile) to identify outlier samples based on number of reads and number of non-zero genes. 2 recommended 
    - Output:
        - Outliers.xlsx: List of samples selected as outliers based on the QC of sequencing depth and number of non-zero genes
        - QC_Data.xlsx: File with count matrix and metadata without selected outliers
        - QC_Prenormalization: QC plots
* PRE_FILTERING: Filtering genes depending on different criteria 
    - Params:
        - Output_file_path: Directory path for storing results
        - PRE_FILTER: PREV (Prevalence of gene across samples is used as a prefiltering criterion)
        - Prev_perc: Threshold used during gene filtering. Genes present in less than Prev_percent samples will be removed
    - Output:
        - Prefilter_Data.xlsx: File with Count matrix with gene prefiltering applied to count matrix
* DESEQ_NORM: DESeq normalization and PCA plots
    - Params:
        - Output_file_path: Directory path for storing results
        - QCNORM: Whether the normalization is made before or after the QC based on gene expression behavior. PRE_QCNORM (Normalization before QC), POST_QCNORM (Normalization after QC)
        - CoarseConditions: Features from metadata to build design matrix of DESeq normalization. These features will be used to select comparisons
        - METHOD_NORM: Normalization method. 'Standard' (Normalization obtained from DESeq object) or 'Standard_Parallel' (Standard DESeq2 method but using parallelization with clustermq) 
        - VST_FILTER: Whether to apply Variance Stabilizing transformation to the normalized matrix. VST_ON or VST_OFF. VST_ON recommended
        - SVD_FILTER: Whether to apply SVD truncation to eliminate noisy PCA direction contributions. SVD_OFF or SVD_ON
        - PlotPCA: Whether to produce PCA plots. 'PCA_PLOT'
        - Control_Neg_PCA: Negative control to show in PCAs. Use "" if only samples associated to a particular treatment are desired in the plot
        - Control_Pos_PCA: Positive control to show in PCA. Use "" if there is none.
        - MEM_MB: Memory parameter for parallelization with register_dopar_cmq from clustermq package. Only used when METHOD_NORM = 'Standard_Parallel'
    - Output:
        - DESeq_Norm.RData: DESeq2 object with normalization and DEG information
        - Norm_Data.xlsx: Normalized count matrix
        - PCA_PLOTS: PCA plots using Norm_Data per treatment and for all samples
          
* QC_POSTNORMALIZATION: Outlier detection based on gene expression and within replicate correlation
    - Params:
        - Output_file_path: Directory path for storing results
        - OUTLIER_FILTER: Method to detect outliers. GENTLE_REP (Method based on correlation within replicates only)
        - RepFeatures_path: Metadata features to identify biological replicates 
   
    - Output:
       - Outlier_DF.xlsx:
            - Bad_Samples_cond: Samples associated to conditions where correlation among replicates is significantly low
            - Outlier_samples: Samples that correlate very poorly with all its replicates
            - Outlier_Selected: Bad_Samples + Outlier_samples. Samples to be removed from the rest of the analysis.
        - Correlation_within_replicates.pdf: Plots showing correlation within replicates across all conditions
        - PCA_WITH_OUTLIERS: PCA plots showing selected outliers
  
* DEG_FUNCTION: Finds differentially expressed genes
    - Params:
        - Output_file_path: Directory path for storing results
        - List_contrasts_Path: Path to .xlsx file containing a list with pairs of conditions for which DEG should be estimated. Columns have to be named "treat" and "untreat"
        - DEG_Method: Method to find DEG. Either 'DESeq' or 'T-Test'
        - MH_Method: Multiple hypothesis testing method. Either 'BH' or 'High_Cr'(high criticism)
        - AlphaHC: Significance threshold for high criticism method. A number between 0 and 1
        - padjval: P adjusted value threshold for significance. A number between 0 and 1. 0.2 is recommended
        - LogFoldThrs_VolPlot: LogFold threshold to use in Volcano plot and to select DEG for pathway enrichment
    - Output: A folder per molecule
        - DEG_MXXX_.xlsx: List containing DEG for every comparison containing the selected molecule
        - Volcano plots for every comparison

* PWAY_ENRICHMENT: Finds enriched pathways for each of the comparisons found in DEG_FUNCTION
       
    - Params:
        - Output_file_path: Directory path for storing results
        - WITH_REVERSE: Whether enrichment using reverse filter should also be estimated. If WITH_REVERSE is TRUE, enriched pathways are calculated using only genes regulated in the healthy direction (REVERSE = TRUE), as well as all DEG (REVERSE = FALSE)
        - REVERSE_DRUG: Name of treatment/state to use as healthy
        - REVERSE_TIME: Collection time in case healthy control is only measured at a particular time
        - padjval: P adjusted value threshold for significance of DEG that will be used in pathway enrichment. A number between 0 and 1. 0.2 is recommended
        - LogFoldThrs: LogFold threshold use for selecting DEGs that will be used in pathway enrichment
        - Pway_qvalThrs:  P adjusted value threshold for significance of enriched pathways. A number between 0 and 1. 0.2 is recommended
        - ORGANISM: Either "Mouse" or "Human"
        - Target_list_path: Path to target list candidates from chemoproteomics and chemoinformatics
        - DirPipeline_Data: Path to directory where KEGG, REACTOME and GO gene sets are stored
        - GeneDescription_path: Path to where the file containing gene description is stored
        - TimeToRemove: When an entire category of controls at a given point is removed, use this parameter to remove all samples at this time point, given that no comparison can be made at this time point if there are no controls
     
* PWAYS_INTEGRATION: Creates summary data frames across compared conditions for pathway enrichment and data frames for over-representation of chemoproteomics candidate targets

    - Params:
         - Output_file_path: Directory path for storing results
         - Target_list_path: Path to target list candidates from chemoproteomics and chemoinformatics
         - ORGANISM: Either "Mouse" or "Human"
         - DirPipeline_Data: Path to directory where KEGG, REACTOME and GO gene sets are stored
         - Control: Name of control used for comparisons
    

          
        
  
    

