library(rmarkdown)
library(tidyverse)
library(readxl)


mydrugs = c("M082")
cell_line = "THP-1"
stimulation = "LPS"
scinamicNum = 'EXP007'
baseDir = "/fsx/home/crivera/Novogene/RESULTS_FROM_PIPELINE/221101_EXP007_BCH003/240708_Results_THP1_LPS"
DirDataForPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline"
DirPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R"
dataDir =  "/fsx/home/crivera/Novogene/RESULTS_FROM_PIPELINE/221101_EXP007_BCH003/240708_Results_THP1_LPS/DESEQ_NORM"
Contrast_path ="/fsx/home/crivera/Novogene/RESULTS_FROM_PIPELINE/221101_EXP007_BCH003/Data_Files_Pipeline/List_contrasts.xlsx"
Chemo_path =  "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/Target_list_240728_w_InvitroTargets.RDS"
Cell_Dict_path = "/fsx/home/crivera/Novogene/RESULTS_FROM_PIPELINE/221101_EXP007_BCH003/Data_Files_Pipeline/CellLine_Dict_Chemo.csv"

Name_folder_path =  file.path(baseDir, "REPORTS")
dir.create(Name_folder_path, recursive = TRUE)

for (x in mydrugs){

    rmarkdown::render(input = paste(DirPipeline,"DEG_ENRICH_perdrug.Rmd", sep = '/')  ,
                      params = list(drug = x,
                                    cell_line = cell_line,
                                    stimulation = stimulation,
                                    scinamicNum = scinamicNum,
                                    baseDir = baseDir,
                                    DirDataForPipeline = DirDataForPipeline,
                                    DirPipeline = DirPipeline,
                                    dataDir = dataDir,
                                    Contrast_path = Contrast_path,
                                    Chemo_path = Chemo_path,
                                    Cell_Dict_path = Cell_Dict_path),clean = TRUE,
                      output_file = file.path(Name_folder_path,paste(format(Sys.time(), '%y-%m-%d'), cell_line,stimulation, x, "DEG_ENRICH", '.html',sep = '_')))

}
