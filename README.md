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

Two reports are generated:
* QC_REPORT:
 
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
       
    
The following functions have been used throughout the pipeline



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
     

    
