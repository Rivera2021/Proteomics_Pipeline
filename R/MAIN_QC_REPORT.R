library(rmarkdown)
library(tidyverse)
library(readxl)


cellline = "K562"
scinamicNum = 'A_2024_0142'
baseDir = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0142/Results_240905"
DirDataForPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline"
DirPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R"
Path_metrics = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0142/Data/MultiQC.tsv"
QCNORM = "PRE_QCNORM"
TPM_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0142/Data/TPM.tsv"
ControlName = "DMSO"
Batch_variable = ""
ExpDescription = "The protocol used for this experiment was Watchmaker"
Prev_perc = 0.2
MEMORY_MB = 4000
Corr_plot_features = c("Cell_line", "Treatment","Treatment_conc_uM", "Treatment_time_hrs", "Stimulant_used","Sample_name")
CoarseCondition =  c("Cell_line", "Treatment","Treatment_conc_uM", "Treatment_time_hrs", "Stimulant_used")
Gene_annotation = ""
Gene_ensemble = TRUE




Name_folder_path =  file.path(baseDir, "REPORTS")
dir.create(Name_folder_path, recursive = TRUE)

# THe correlation plots are not being plotted using this function. We need to solve that using the childs or something else.

rmarkdown::render(input = paste(DirPipeline,"QC_REPORT.Rmd", sep = '/')  ,
                      params = list(cellline =  cellline,
                                    scinamicNum = scinamicNum,
                                    baseDir =  baseDir,
                                    DirDataForPipeline = DirDataForPipeline,
                                    DirPipeline = DirPipeline,
                                    Path_metrics = Path_metrics,
                                    QCNORM =  QCNORM,
                                    TPM_path = TPM_path,
                                    ControlName = ControlName,
                                    Batch_variable = Batch_variable,
                                    ExpDescription = Batch_variable,
                                    Prev_perc =  Prev_perc,
                                    MEMORY_MB = MEMORY_MB,
                                    Corr_plot_features = Corr_plot_features,
                                    CoarseCondition =  CoarseCondition,
                                    Gene_annotation = Gene_annotation,
                                    Gene_ensemble = Gene_ensemble),clean = TRUE,
                      output_file = file.path(Name_folder_path,paste(format(Sys.time(), '%y-%m-%d'), cellline, "QC_REPORT", '.html',sep = '_')))


