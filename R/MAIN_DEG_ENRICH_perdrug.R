library(rmarkdown)
library(tidyverse)
library(readxl)


cell_line = "K562"
stimulation = "TNFa"
scinamicNum = 'A-2024-0085'
baseDir = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0085/10_million/Results_240818"
DirDataForPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline"
DirPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R"
dataDir =  "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0085/10_million/Results_240818/DESEQ_NORM"
Contrast_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0139_AND_A_2024_0140/Data/List_contrasts_HBEC5i_Time_4hrs.xlsx"
Chemo_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/Target_list_240814_w_InvitroTargets.RDS"
Cell_Dict_path = ""
AdjustDeSeq = ""


Name_folder_path =  file.path(baseDir, "REPORTS")
dir.create(Name_folder_path, recursive = TRUE)

Metadata <-read_xlsx(file.path(dataDir, "Norm_Data.xlsx"), sheet = 'Metadata')
drugs = unique(Metadata$Treatment)
drugs = drugs[drugs!= 'DMSO']
drugs = drugs[drugs!= "NONE"]


for(x in drugs){
    rmarkdown::render(input = paste(DirPipeline,"DEG_ENRICH_perdrug.Rmd", sep = '/')  ,
                      params = list(drug = x,
                                    cell_line = cell_line,
                                    scinamicNum = scinamicNum,
                                    baseDir = baseDir,
                                    DirDataForPipeline = DirDataForPipeline,
                                    DirPipeline = DirPipeline,
                                    dataDir = dataDir,
                                    Contrast_path = Contrast_path,
                                    Chemo_path = Chemo_path,
                                    AdjustDeSeq = AdjustDeSeq),clean = TRUE,
                      output_file = file.path(Name_folder_path,paste(format(Sys.time(), '%y%m%d'), cell_line,stimulation, x, "DEG_ENRICH", '.html',sep = '_')))

}

