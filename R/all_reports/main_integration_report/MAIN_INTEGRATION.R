library(rmarkdown)
library(tidyverse)
library(readxl)




ConditionsComparisons_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0142_AND_A_20204_0143/Data/Integration_sheet_50_Percent_rapa.xlsx"
CoarseCond =  "K562_NONE_rapamycin"
baseDir = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0142_AND_A_20204_0143/Results_240826"
DirPipeline =  "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R"
DirDataForPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline"
Chemo_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/Target_list_c20_240826.RDS"
scinamicNums = "A-2024-0142"

# ScinamicNumbers = unique(ConditionsComparisons$`Experiment ID`)
# ConComp_tabl <-read_xlsx(ConditionsComparisons_path)
# ConComp_tabl$CoarseCond = paste(ConComp_tabl$Cell_line, ConComp_tabl$Stimulation, ConComp_tabl$Treatment, sep = '_')
# CoarseCond = unique(ConComp_tabl$CoarseCond)


Name_folder_path =  file.path(baseDir, "REPORTS")
dir.create(Name_folder_path, recursive = TRUE)


rmarkdown::render(input = paste(DirPipeline,"INTEGRATION_REPORT.Rmd", sep = '/')  ,
                      params = list(ConditionsComparisons_path = ConditionsComparisons_path,
                                    CoarseCond = CoarseCond,
                                    scinamicNums = scinamicNums,
                                    baseDir = baseDir,
                                    DirDataForPipeline = DirDataForPipeline,
                                    Chemo_path = Chemo_path,
                                    DirPipeline = DirPipeline),clean = TRUE,
                      output_file = file.path(Name_folder_path,paste(format(Sys.time(), '%y%m%d'),CoarseCond, "INTEGRATION_REPORT", '.html',sep = '_')))


