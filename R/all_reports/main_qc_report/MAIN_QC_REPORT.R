library(rmarkdown)
library(tidyverse)
library(readxl)

#source("/fsx/home/john/projects/empress/old_2022_runs/compound2/230404_EXP008_BCH003/params_qc.txt")

# Create contrasts sheets for all the experiments in folder ~/MULTI_OMICS/ChemoBetter_results/Proteomics/2025

#baseDir_all <- "/fsx/home/crivera/MULTI_OMICS/ChemoBetter_results/Proteomics/2025_Corrected/omics-proteomics"
#baseDir_raw = "/fsx/home/crivera/MULTI_OMICS/ChemoBetter_results/Proteomics/2025_Raw_Corrected"
baseDir_all <- "/fsx/home/crivera/MULTI_OMICS/ChemoBetter_results/Proteomics/2025_Corrected_exp155B"
baseDir_raw = "/fsx/home/crivera/MULTI_OMICS/ChemoBetter_results/Proteomics/2025_Raw_Corrected_exp155B"
Exps = list.files(baseDir_all)
OutputFolder = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/PROTEOMICS/ChemoBetter/Results_252224_PROTEOMICS_exp155B"
DirPipeline = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/all_reports/main_qc_report"
Alias_path =  "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/mol_ids_chemobetter.csv"
Uniprot_map_db = "/fsx/home/crivera/Data_Bases/UNIPROT_GeneInfo/hgnc_uniprot_mapping.txt"
workst = "Proteomics"

for(exp in Exps){
#elnid: uncomment with new format of organization
x = strsplit(exp, "-")[[1]]
exp_real = paste(x[1:(length(x)-1)], collapse = '-')

dir.create(file.path(OutputFolder, exp, "REPORTS"),  recursive = TRUE)
Rawfiles = list.files(file.path(baseDir_raw, paste("elnid=", exp_real,sep = ""), paste("expid=",exp,sep = "")))
Extract_Metadata_name = Rawfiles[grep("Extract Proteomics", Rawfiles)]
protdisc_name = Rawfiles[grep("Proteins.txt", Rawfiles)]
PSM_file_name = Rawfiles[grep("PSMs.txt", Rawfiles)]



rmarkdown::render(input = paste(DirPipeline,"QC_REPORT.Rmd", sep = '/')  ,
                  params = list(scinamicNum = exp,
                                baseDir =  file.path(OutputFolder, exp),
                                DirPipeline = DirPipeline,
                                Alias_path = Alias_path,
                                Uniprot_map_db = Uniprot_map_db,
                                Metadata_gr_path = file.path(baseDir_all, exp, paste(exp, "contrast_metadata.csv", sep = '_')),
                                Contrast_def_path = file.path(baseDir_all, exp, paste(exp, "contrast_definitions.csv", sep = '_')),
                                workstream = workst,
                                Extract_Metadata_path = file.path(baseDir_raw, paste("elnid=", exp_real ,sep = ""),paste("expid=",exp,sep = ""), Extract_Metadata_name),
                                Proteins_protdisc_path = file.path(baseDir_raw, paste("elnid=", exp_real ,sep = ""),paste("expid=",exp,sep = ""), protdisc_name),
                                Contrast_groups_path = file.path(baseDir_all, exp, paste(exp, "contrast_groups.csv", sep = '_')),
                                PSM_path = file.path(baseDir_raw, paste("elnid=", exp_real ,sep = ""),paste("expid=",exp,sep = ""), PSM_file_name)
                  ),
                  clean = TRUE,
                  output_file = file.path(OutputFolder, exp, "REPORTS",paste(format(Sys.time(), '%y-%m-%d'),exp, "QC_REPORT", '.html',sep = '_')))

}


