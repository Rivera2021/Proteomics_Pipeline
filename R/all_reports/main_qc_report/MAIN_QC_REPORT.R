library(rmarkdown)
library(tidyverse)
library(readxl)

cellline = "In vivo"
scinamicNum = 'EXP000033'
baseDir = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Results_241115"
DirDataForPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline"
DirPipeline =  "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/all_reports/main_qc_report"
Metrics_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Data/MultiQC.tsv"
Metadata_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Data/241030_EXP000033_metadata_v3_MB.xlsx"
CountM_path =  "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Data/salmon.merged.gene_counts.tsv"
SampleName_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Data/SampleName.csv"
Experimental_design_path = ""
Column_match = "ID"
Organism_type = "In_vivo" # Either In_vivo or In_vitro
QCNORM: "PRE_QCNORM"
TPM_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Data/salmon.merged.gene_tpm.tsv"
ControlName = "Vehicle-DMSO"
Batch_variable= ""
ExpDescription = "The protocol used for this experiment was Watchmaker"
Prev_perc = 0.2
MEMORY_MB = 4000
Corr_plot_features = c("Treatment","Treatment_conc", "Treatment_time_hrs", "Stimulant_used","Sample_name")
CoarseCondition =  c("Treatment","Treatment_conc", "Treatment_time_hrs", "Stimulant_used")
Alias_path = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/mol_ids.csv"
PRE_FILTER = "PREV"



Name_folder_path =  file.path(baseDir, "REPORTS")
dir.create(Name_folder_path, recursive = TRUE)

# THe correlation plots are not being plotted using this function. We need to solve that using the childs or something else.

rmarkdown::render(input = paste(DirPipeline,"QC_REPORT.Rmd", sep = '/')  ,
                      params = list(cellline =  cellline,
                                    scinamicNum = scinamicNum,
                                    baseDir =  baseDir,
                                    DirDataForPipeline = DirDataForPipeline,
                                    DirPipeline = DirPipeline,
                                    Metrics_path = Metrics_path,
                                    Metadata_path = Metadata_path,
                                    CountM_path =CountM_path,
                                    SampleName_path = SampleName_path,
                                    Experimental_design_path = Experimental_design_path,
                                    Column_match = Column_match,
                                    Organism_type = Organism_type,
                                    QCNORM =  QCNORM,
                                    TPM_path = TPM_path,
                                    ControlName = ControlName,
                                    Batch_variable = Batch_variable,
                                    ExpDescription = Batch_variable,
                                    Prev_perc =  Prev_perc,
                                    MEMORY_MB = MEMORY_MB,
                                    Corr_plot_features = Corr_plot_features,
                                    CoarseCondition =  CoarseCondition,
                                    Alias_path = Alias_path,
                                    PRE_FILTER = PRE_FILTER),clean = TRUE,
                      output_file = file.path(Name_folder_path,paste(format(Sys.time(), '%y-%m-%d'), cellline, "QC_REPORT", '.html',sep = '_')))


