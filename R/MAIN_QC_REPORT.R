library(rmarkdown)
library(tidyverse)
library(readxl)



cellline =  "K562, HBEC-5i"
scinamicNum = 'A_2024_0139, A_2024_0140, A_2024_0142, A_2024_0143 '
baseDir = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0139_AND_A_2024_0140_A_2024_0142_A_2024_0143/Results_240809"
DirDataForPipeline ="/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline"
DirPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R"
CellLine_Dict_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0139_AND_A_2024_0140_A_2024_0142_A_2024_0143/Data/CellLine_Dict.csv"
Path_metrics = ""
QCNORM = "PRE_QCNORM"
TPM_path = " "
ControlName = "DMSO"
Batch_variable = ""
ExpDescription = "The protocol used for this experiment was watchmaker"
Prev_perc = 0.2
MEMORY_MB = 4000







Name_folder_path =  file.path(baseDir, "REPORTS")
dir.create(Name_folder_path, recursive = TRUE)

# THe correlation plots are not being plotted using this function. We need to solve that using the childs or something else.

rmarkdown::render(input = paste(DirPipeline,"QC_REPORT.Rmd", sep = '/')  ,
                      params = list(cellline = cellline,
                                    scinamicNum = scinamicNum,
                                    baseDir = baseDir,
                                    DirDataForPipeline = DirDataForPipeline,
                                    DirPipeline = DirPipeline,
                                    CellLine_Dict_path = CellLine_Dict_path,
                                    Path_metrics = Path_metrics,
                                    QCNORM = QCNORM,
                                    TPM_path = TPM_path,
                                    ControlName = ControlName,
                                    Batch_variable = Batch_variable,
                                    ExpDescription = ExpDescription,
                                    Prev_perc = Prev_perc,
                                    MEMORY_MB = MEMORY_MB),clean = TRUE,
                      output_file = file.path(Name_folder_path,paste(format(Sys.time(), '%y-%m-%d'), cellline, "QC_REPORT", '.html',sep = '_')))


