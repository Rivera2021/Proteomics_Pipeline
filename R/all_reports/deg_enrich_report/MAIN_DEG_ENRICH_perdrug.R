library(rmarkdown)
library(tidyverse)
library(readxl)


cell_line = "In_vivo"
stimulation = "TNBS"
scinamicNum = 'EXP000033'
baseDir = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Results_241115"
DirDataForPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline"
DirPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/all_reports/deg_enrich_report"
dataDir =  "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Results_241115/DESEQ_NORM"
Contrast_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Data/List_contrasts_Invivo_v2.xlsx"
Chemo_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/Target_list_241112_w_InvitroTargets_knInh.RDS"
expr_gene_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Results_241115/PRE_FILTERING/Expressed_genes.xlsx"
Experimental_design_path = ""
Cell_Dict_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Data/CellLine_Dict_Chemo.csv"
AdjustDeSeq = ""
Organism = "Mouse"
padj_thr_gene = 0.2
logFC_thrs_gene =  0.5


Name_folder_path =  file.path(baseDir, "REPORTS")
dir.create(Name_folder_path, recursive = TRUE)

Metadata <-read_xlsx(file.path(dataDir, "Norm_Data.xlsx"), sheet = 'Metadata')
drugs = unique(Metadata$Treatment)
drugs = drugs[drugs!= 'Vehicle-DMSO']
drugs = drugs[drugs!= "NONE"]
drugs = "M255"

for(x in drugs){
    rmarkdown::render(input = paste(DirPipeline,"DEG_ENRICH_perdrug.Rmd", sep = '/')  ,
                      params = list(drug = x,
                                    stimulation = stimulation,
                                    cell_line = cell_line,
                                    scinamicNum = scinamicNum,
                                    baseDir = baseDir,
                                    DirDataForPipeline = DirDataForPipeline,
                                    DirPipeline = DirPipeline,
                                    dataDir = dataDir,
                                    Contrast_path = Contrast_path,
                                    Chemo_path = Chemo_path,
                                    Experimental_design_path = Experimental_design_path,
                                    Cell_Dict_path = Cell_Dict_path,
                                    AdjustDeSeq = AdjustDeSeq,
                                    Organism =  Organism,
                                    padj_thr_gene = padj_thr_gene,
                                    logFC_thrs_gene = logFC_thrs_gene),clean = TRUE,
                      output_file = file.path(Name_folder_path,paste(format(Sys.time(), '%y%m%d'), cell_line,stimulation, x, "DEG_ENRICH", '.html',sep = '_')))

}
