library(rmarkdown)
library(tidyverse)
library(readxl)


cellline = "THP-1"
scinamicNum = 'EXP004'
baseDir = "/fsx/home/crivera/Novogene/RESULTS_FROM_PIPELINE/221101_EXP007_BCH003/240708_Results_THP1_LPS"
DirDataForPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline"
DirPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R"
CellLine_Dict_path = "/fsx/home/crivera/Novogene/RESULTS_FROM_PIPELINE/221101_EXP007_BCH003/Data_Files_Pipeline/CellLine_Dict.csv"
Path_metrics = "/fsx/home/crivera/Novogene/RESULTS_FROM_PIPELINE/221101_EXP007_BCH003/Data_Files_Pipeline/multiqc_general_stats.txt"
QCNORM = "PRE_QCNORM"
TPM_path = "3prime"
ControlName = "DMSO"
ExpDescription = "The cells of interest were seeded in the 96 well plates with planned layout. Four hours before adding the stimulation, cells were pretreated with either MOLs or DMSO. After stimulation for another 4 hours, cells were lysated in TCL or RLT buffer. Cell lysis can be store to -80C.\nRNA was purified by adding 2.2 volume of RNA cleanup beads. DNase digestion was used to remove genomic DNA contamination. 20 ng RNA was used for reverse transcription using mRNA specific oligo dT primers with barcode and UMI. The library was prepared using Nextera XT DNA Library Preparation Kit according to the instruction. NovaSeq6000 seq platform was used to sequence the pooled libraries. Pair-end sequence and 25-8-0-151 sequence cycle were used. The raw data were trimmed to 16-8-0-150 and pre-demulitplexed by Novogene."

Name_folder_path =  file.path(baseDir, "REPORTS")
dir.create(Name_folder_path, recursive = TRUE)



    rmarkdown::render(input = paste(DirPipeline,"QC_REPORT.Rmd", sep = '/')  ,
                      params = list(cellline = cellline,
                                    scinamicNum = scinamicNum,
                                    baseDir = baseDir,
                                    DirDataForPipeline = DirDataForPipeline,
                                    DirPipeline = DirPipeline,
                                    CellLine_Dict_path = CellLine_Dict_path,
                                    Path_metrics = Path_metrics,
                                    QCNORM = QCNORM,
                                    Cell_Dict_path = Cell_Dict_path),clean = TRUE,
                      output_file = file.path(Name_folder_path,paste(format(Sys.time(), '%y-%m-%d'), cell_line,stimulation, x, "DEG_ENRICH", '.html',sep = '_')))



