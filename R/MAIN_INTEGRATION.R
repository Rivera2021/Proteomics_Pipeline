library(rmarkdown)
library(tidyverse)
library(readxl)



ConditionsComparisons_path = "/fsx/home/crivera/Novogene/RESULTS_FROM_PIPELINE/221101_EXP007_BCH003/Data_Files_Pipeline/Integration_Comparisons.xlsx"
baseDir = "/fsx/home/crivera/Novogene/RESULTS_FROM_PIPELINE/221101_EXP007_BCH003/240708_Results_THP1_LPS"
DirPipeline =  "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R"
DirDataForPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline"

ScinamicNumbers = unique(ConditionsComparisons$`Experiment ID`)
ConComp_tabl <-read_xlsx(ConditionsComparisons_path)
ConComp_tabl$CoarseCond = paste(ConComp_tabl$Cell_line, ConComp_tabl$Stimulation, ConComp_tabl$Treatment, sep = '_')
CoarseCond = unique(ConComp_tabl$CoarseCond)


Name_folder_path =  file.path(baseDir, "REPORTS")
dir.create(Name_folder_path, recursive = TRUE)


rmarkdown::render(input = paste(DirPipeline,"INTEGRATION_REPORT.Rmd", sep = '/')  ,
                      params = list(ConditionsComparisons_path = ConditionsComparisons_path,
                                    CoarseCond = CoarseCond,
                                    scinamicNums = ScinamicNumbers,
                                    baseDir = baseDir,
                                    DirDataForPipeline = DirDataForPipeline,
                                    DirPipeline = DirPipeline),clean = TRUE,
                      output_file = file.path(Name_folder_path,paste("INTEGRATION_REPORT", '.html',sep = '_')))



