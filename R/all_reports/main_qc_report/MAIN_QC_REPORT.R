library(rmarkdown)
library(tidyverse)
library(readxl)

# experimental params
scinamicNum = 'EXP000064'
cellline = "K562.1 & THP1.1"
Experimental_design_path = ""
Column_match = "ID"
Organism_type = "In_vitro" # Either In_vivo or In_vitro
ControlName = "Vehicle-DMSO"
Batch_variable= ""
ExpDescription = "The protocol used for this experiment was Watchmaker"

# data and metadata locations
expDir = "/fsx/home/john/projects/empress/EXP000064"
baseDir = file.path(expDir, "Results_250301/")
Metrics_path = file.path(expDir, "rawdata/multiqc/star_salmon/multiqc_report_data/")
Metadata_path = file.path(expDir, "meta/EXP000064_metadata.xlsx")
CountM_path =  file.path(expDir, "rawdata/star_salmon/salmon.merged.gene_counts.tsv")
SampleName_path = file.path(expDir, "meta/SampleName.csv")
TPM_path = file.path(expDir, "rawdata/star_salmon/salmon.merged.gene_tpm.tsv")

# report params
Prev_perc = 0.2
MEMORY_MB = 4000
Corr_plot_features = c("Treatment","Treatment_conc", "Treatment_time_hrs", "Stimulant_used","Sample_name")
CoarseCondition =  c("Cell_line", "Treatment","Treatment_conc", "Treatment_time_hrs", "Stimulant_used")
PRE_FILTER = "PREV"
QCNORM = "PRE_QCNORM"

# pipeline specific params
DirDataForPipeline = "/fsx/home/john/repos/Transcriptomics_Pipeline/data/Data_for_pipeline"
DirPipeline =  "/fsx/home/john/repos/Transcriptomics_Pipeline/R/all_reports/main_qc_report"
Alias_path = "/fsx/home/john/repos/Transcriptomics_Pipeline/data/Data_for_pipeline/mol_ids.csv"



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
                                PRE_FILTER = PRE_FILTER),
                  clean = TRUE,
                  output_file = file.path(Name_folder_path,paste(format(Sys.time(), '%y-%m-%d'), cellline, "QC_REPORT", '.html',sep = '_')))


