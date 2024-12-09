# Transcriptomics_Pipeline

This project is used to create the Nextflow pipeline for bulk RNAseq transcriptomics data analysis. 

Started on 2023-10-27.

## Directory structure

* **R** - Resuable R code (functions etc.)
* **data** - Raw data used for analysis

## Data

* Data_for_pipeline: Contains different type of data that is used along the pipeline
* In_vivo: Contains unit tests data sets, raw count matrix and metadata from an in vivo study. It also contains a few files needed in the pipeline.



## Code

Three reports are generated:
* QC_REPORT: Contains technical QC, PCA plots, and correlations plots. Uses the following files
    - MAIN_QC_REPORT.R: Renders QC_REPORT.Rmd 
        - This is the list of parameters required to generate this report:
            - cellline:  Parameter used ONLY in the title of the QC reports to identify the cell lines included in the QC analysis.
            - scinamicNum: Either scinamic number or Experiment number. Will be ONLY used in the title of the QC report
            - baseDir: Directory where all the files will be stored
            - DirDataForPipeline: Directory where pipeline files are being stored
            - DirPipeline: Path where the R scripts and .rmd files for QC are being stored
            - Metrics_path: Path to af file named “multiqc_general_stats.txt” that nf-core rnaseq pipeline generates within the multiqc folder.
            - Metadata_path: Path to the metadata excel file
            - CountM_path: Path to the Count matrix. This is the file named as “salmon.merged.gene_counts.tsv” that nf-core rnaseq pipeline within the star_salmon folder. Should have a column with gene_ID and another column with gene_name (using symbol notation). The column names in this file should match one of the columns from the Metadata file
            - SampleName_path: Path to a file that contains the names in the metadata that should be used to name the samples for the rest of the analysis. Should end up in a unique combination for every sample. Possible options are: Cell_line, Treatment, Treatment_conc,  Stimulant_used, Treatment_time_hrs, Replicate.
            - Experimental_design_path: Path to a png image with the experimental design
            - Column_match: Column from metadata to match names in count matrix.
            - Organism_type: Either “In_vivo” or “In_vitro”
            - QCNORM: Either “PRE_QCNORM” or  “POST_QCNORM” (use this when want to avoid including the customized list of outliers for the rest of the analysis)
            - ControlName: Name of the Control. For instance: “DMSO”, “Vehicle”
            - Batch_variable: single name from the metadata that wants to be batch corrected. For example: Plate
            - ExpDescription: A text descrbing the experiment
            - PRE_FILTER: Prefiltering type use to get rid of low express genes. Either “PREV” based on prevalence of the gene across samples, or “LOW_EXPR” getting rid of genes that have less number of samples than the smallest group size with at least 10 reads. 
            - Prev_perc: Percent prevalence of a genes across all samples. Can be selected based on the percent of the smallest group size. 
            - MEMORY_MB: Memory usage for using the parallelize version of DeSeq using slurm
            - Corr_plot_features: Features from the metadata to be added in the correlation plots. Possible values are: Cell_line, Treatment, Treatment_conc,  Stimulant_used, Treatment_time_hrs, RQN.
            - CoarseCondition: Features from the metadata to select groups of samples for contrasts. Should uniquely categorize samples within conditions to be compared. Possible values are: "Treatment","Treatment_conc", "Treatment_time_hrs"
            - Alias_path: Path mapping compound names from scinamic to its alias.
            - TPM_path: Path to the TPM matrix. This is the file named as “salmon.merged.gene_tpm.tsv” that nf-core rnaseq pipeline outputs within the star_salmon folder
            - outliers_path: Path to the outliers that want to be removed if used QCNORM = “POST_QCNORM”. If QCNORM = “PRE_QCNORM” leave as "". This is a file that contains the user-selected outliers after looking at the QC_REPORT generated with all samples, QCNORM = “PRE_QCNORM”. It is a .xlsx file. Use the example generated from the file Outliers.xlsx within QC_PRENORMALIZATION folder, or the file within In_vitro_EXP72. Notice that outlier samples should be identified using the Sample_name column from the metadata generated in QC_PRENORMALIZATION.
            

        - The following functions are used:
            - IMPORT_DATA
            - QC_PRENORMALIZATION_CL_V2
            - PRE_FILTERING
            - DESEQ_NORM
            
 
* DEG_ENRICH: For every cell line and every molecule a report is generated. It contains DEGs and GSEA enriched pathways. Uses the following files
    - MAIN_DEG_ENRICH_perdrug.R: Renders MAIN_DEG_ENRICH_perdrug.R
    
        - This is the list of parameters required to generate this report:
        
            - Drug:  Treatment for which report is made. Should use the alias name
            - cell_line: Parameter used ONLY in the title of the DEG_ENRICH reports to identify the cell lines included in the DEG_ENRICH analysis.
            - Stimulation:  Parameter used ONLY in the title of the DEG_ENRICH reports to identify the stimulations described in this report.
            -  scinamicNum: Either scinamic number or Experiment number. Will be ONLY used in the title of the report.
            -  baseDir: Directory where all the files will be stored
            -  DirDataForPipeline: Directory where pipeline files are being stored
            -  DirPipeline: Path where the R scripts and .rmd files for DEG and enrichment are being stored
            -  dataDir: Directory where the deseq object with normalized matrix and filter count matrix can be found. Should be within the DESE_NORM folder generated in the QC analysis
            - Contrast_path: Path to excel file determining the contrasts desired. Only those related to the molecule associated to this report or controls contrats swill be display in this report. File should have as column names ‘treat’ and ‘untreat’. The analysis will be made treat_vs_untreat. The way to identify conditions should be using the coarse_condition names generated in the QC analysis.
            - Chemo_path: Path to excel file containing chemoproteomics information
            - expr_gene_path: Path to file named “Expressed_genes.xlsx ” generated in the QC analysis (PRE_FILTERING folder). This file identifies genes that are expressed in the cell line.
            - Experimental_design_path: Path to a png image with the experimental design
            - Cell_Dict_path: path to .csv file mapping cell lines in transcriptomics metadata and cell lines in chemoproteomics metadata
            - AdjustDeSeq: single name from the metadata that wants to be batch corrected. For example: Plate
            - Organism: Either “Mouse” or “Human”
            - padj_thr_gene: Threshold in adjusted pvalue to identify differentially expressed genes
            - logFC_thrs_gene: Threshold in Log fold-change to identify differentially expressed genes

    - The following functions are used:
        - DEG_FUNCTION_DA.R
    
    
* MAIN_INTEGRATION (In progress): Contains integration with chemoproteomics across selected conditions. Uses the following files
    - MAIN_INTEGRATION.R: Renders INTEGRATION_REPORT.Rmd


The following functions have been used throughout the pipeline


* IMPORT_DATA.R (Not used until protocol is established): Validates consistency between count matrix and Metadata
    - Params:
        - Metadata_path: Path to Metadata in .xlsx format
        - Count_path: Path to matrix with read count in .tsv format
        - Match_feature: Column name from metadata that matched the sample names used in the count matrix
        - Output_file_path: Directory path for storing results
        - SampleNamepath: Path to .csv file containing the metadata column names that will be used to uniquely identify a sample using biological features
    - Output:
        - Import_Data.xlsx: File with count matrix and metadata after a few modifications needed to start the pipeline
* QC_PRENORMALIZATION_CL_V2: QC based on the quality of the sequencing
    - Params:
        - Output_file_path: Directory path for storing results
        - nMust: Number of IQRs below Q1(25% quantile) to identify outlier samples based on number of reads and number of non-zero genes. 2 recommended 
        - TPM_path: Path where the TPM file is stored. TPM file is usually output by the nfcore nextflow pipeline
    - Output:
        - Outliers.xlsx: List of samples selected as outliers based on the QC of sequencing depth and number of non-zero genes
        - QC_Data.xlsx: File with count matrix and metadata without selected outliers
        - QC_Prenorm: QC plots
        - TPM_count_expr.RDS: Saves TPM file with average expression across samples
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
        - QCNORM: Whether the normalization is made before or after the QC and outlier detection. PRE_QCNORM (Normalization before QC), POST_QCNORM (Normalization after QC)
        - outliers_path: Path to .xlsx file with outliers names
        - CoarseConditions: Features from metadata to build design matrix of DESeq normalization. These features will be used to select comparisons
        - METHOD_NORM: Normalization method. "Standard"" (Normalization obtained from DESeq object), 'Standard_Parallel' (Parallel mode to use DESeq )
        - VST_FILTER: Whether to apply Variance Stabilizing transformation to the normalized matrix. VST_ON or VST_OFF. VST_ON recommended
        - SVD_FILTER: Whether to apply SVD truncation to eliminate noisy PCA direction contributions. SVD_OFF or SVD_ON
        - PlotPCA (Not used anymore): Whether to produce PCA plots. 'PCA_PLOT'.
        - Control_Neg_PCA (Not used anymore): Negative control to show in PCAs. Use "" if only samples associated to a particular treatment are desired in the plot
        - Control_Pos_PCA (Not used anymore): Positive control to show in PCA. Use "" if there is none.
        - MEM_MB: Memory parameter for parallelization with register_dopar_cmq from clustermq package. Only used when METHOD_NORM = 'Standard_Parallel'

    - Output:
        - DESeq_Norm.RData: DESeq2 object with normalization and DEG information
        - Norm_Data.xlsx: Normalized count matrix
        
          
* QC_POSTNORMALIZATION (Currently not used, but should be included in the QC_REPORT): Outlier detection based on gene expression and within replicate correlation
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
  
* DEG_FUNCTION_DA: Finds differentially expressed genes
    - Params:
        - Output_file_path: Directory path for storing results
        - List_contrasts_Path: Path to .xlsx file containing a list with pairs of conditions for which DEG should be estimated. Columns have to be named "treat" and "untreat"
        - DEG_Method: Method to find DEG. 'DESeq_Cons' (DESeq used only on the samples needed for the commparison)
        - MH_Method: Multiple hypothesis testing method. Either 'BH' or 'High_Cr'(high criticism)
        -  AdjustDeSeq: A character from the metadata features to correct for batch in the DESeq formula

OLD_FUNCTIONS: (Not used in the current pipeline)        
        
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
    

          
        
  
    
