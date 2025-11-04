library(rmarkdown)
library(tidyverse)
library(readxl)

#source("/fsx/home/john/projects/empress/old_2022_runs/compound2/230404_EXP008_BCH003/params_qc.txt")

# Create contrasts sheets for all the experiments in folder ~/MULTI_OMICS/ChemoBetter_results/Proteomics/2025

baseDir_all <- "/fsx/home/crivera/MULTI_OMICS/ChemoBetter_results/Proteomics/2025"
Exps = list.files(baseDir_all)
OutputFolder = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/PROTEOMICS/ChemoBetter"
DirPipeline = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/all_reports/main_qc_report"
Alias_path =  "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/mol_ids.csv"
Uniprot_map_db = "/fsx/home/crivera/Data_Bases/UNIPROT_GeneInfo/hgnc_uniprot_mapping.txt"


for(exp in Exps){

dir.create(file.path(OutputFolder, exp, "REPORTS"),  recursive = TRUE)

rmarkdown::render(input = paste(DirPipeline,"QC_REPORT.Rmd", sep = '/')  ,
                  params = list(scinamicNum = exp,
                                baseDir =  file.path(OutputFolder, exp),
                                DirPipeline = DirPipeline,
                                Alias_path = Alias_path,
                                Uniprot_map_db = Uniprot_map_db,
                                Metadata_gr_path = file.path(baseDir_all, exp, paste(exp, "contrast_metadata.csv", sep = '_')),
                                Contrast_def_path = file.path(baseDir_all, exp, paste(exp, "contrast_definitions.csv", sep = '_'))
                                ),
                  clean = TRUE,
                  output_file = file.path(OutputFolder, exp, "REPORTS",paste(format(Sys.time(), '%y-%m-%d'),exp, "QC_REPORT", '.html',sep = '_')))

}
