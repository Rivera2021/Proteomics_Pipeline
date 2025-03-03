library(rmarkdown)
library(tidyverse)
library(readxl)



scinamicNum = 'EXP000085'
baseDir = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Organoids/EXP000085/Results_250225_PerDonorReport"
DirDataForPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline"
DirPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/all_reports/deg_enrich_report"
dataDir = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Organoids/EXP000085/Results_250225_PerDonorReport/DESEQ_NORM"
Contrast_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Organoids/EXP000085/Data/List_contrasts_perDonor.xlsx"
Chemo_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/Target_list_241210_w_InvitroTargets_knInh.RDS"
expr_gene_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Organoids/EXP000085/Results_250225_PerDonorReport/PRE_FILTERING/Expressed_genes.xlsx"
Experimental_design_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Organoids/EXP000085/Data/Organoids_exp_design_exp85_exp88.png"
Cell_Dict_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Organoids/EXP000085/Data/CellLine_Dict_Chemo.csv"
AdjustDeSeq = ""
Organism = "Human"
padj_thr_gene = 0.2
logFC_thrs_gene = 0
METHOD_NORM = 'Standard_Parallel'
DEG_Method = 'DESeq_Cons'

List_contrasts = read_excel("/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Organoids/EXP000085/Data/List_contrasts_perDonor.xlsx")
Name_folder_path =  file.path(baseDir, "REPORTS")
dir.create(Name_folder_path, recursive = TRUE)

# Metadata <-read_xlsx(file.path(dataDir, "Norm_Data.xlsx"), sheet = 'Metadata')
# drugs = unique(Metadata$Treatment)
# drugs = drugs[!drugs %in% c('DMSO')]
# #drugs = c("dBET6", "TAK-279","Selisistat" )

# Cycle through reportsnum

start_time <- Sys.time()

for(reportnum in unique(List_contrasts$reportNum)){
  # params$drug, params$stimulation and params$cell_line are only used in the title. Therefore they can be concatenation of all the molecules involved
  # in a given report. To find the drug associated to a given contrast inside the report we use the value within the List_contrast_file not the params$drug parameter
    List_contrasts_report = List_contrasts %>% dplyr::filter(reportNum == reportnum)
    cell_line = paste(unique(List_contrasts_report$cell_line), collapse = "_")
    drug= paste(unique(List_contrasts_report$drug), collapse = "_")
    rmarkdown::render(input = paste(DirPipeline,"DEG_ENRICH_perdrug.Rmd", sep = '/')  ,
                      params = list(drug = paste(unique(List_contrasts_report$drug), collapse = "_"),
                                    stimulation = paste(unique(List_contrasts_report$stimulation), collapse = "_"),
                                    cell_line = paste(unique(List_contrasts_report$cell_line), collapse = "_"),
                                    scinamicNum = scinamicNum,
                                    baseDir = baseDir,
                                    DirDataForPipeline = DirDataForPipeline,
                                    DirPipeline = DirPipeline,
                                    dataDir = dataDir,
                                    Contrast_path = Contrast_path,
                                    Chemo_path = Chemo_path,
                                    expr_gene_path = expr_gene_path,
                                    Experimental_design_path = Experimental_design_path,
                                    Cell_Dict_path = Cell_Dict_path,
                                    AdjustDeSeq = AdjustDeSeq,
                                    Organism =  Organism,
                                    padj_thr_gene = padj_thr_gene,
                                    logFC_thrs_gene = logFC_thrs_gene,
                                    ReportNum = reportnum,
                                    METHOD_NORM  = METHOD_NORM ,
                                    DEG_Method = DEG_Method),clean = TRUE,
                      output_file = file.path(Name_folder_path,paste(format(Sys.time(), '%y%m%d'), cell_line, drug,reportnum, "DEG_ENRICH", '.html',sep = '_')))


}

end_time <- Sys.time()

execution_time <- end_time - start_time



