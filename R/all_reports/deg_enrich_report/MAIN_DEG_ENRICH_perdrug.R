library(rmarkdown)
library(tidyverse)
library(readxl)

# Runs enrichments for all experiments in "/fsx/home/crivera/MULTI_OMICS/ChemoBetter_results/Proteomics/2025" at the same time

baseDir_all <- "/fsx/home/crivera/MULTI_OMICS/ChemoBetter_results/Proteomics/2025_Corrected/omics-proteomics"
Exps = list.files(baseDir_all)
OutputFolder = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/PROTEOMICS/ChemoBetter/Results_251121_PROTEOMICS"
DirPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/all_reports/deg_enrich_report"
DirDataForPipeline = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline"
Chemo_path =  "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/Target_list_241210_w_InvitroTargets_knInh.RDS"
Cell_Dict_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/PROTEOMICS/ChemoBetter/Data/CellLine_Dict_Chemo.csv"
Organism = "Human"
padj_thr_gene = 0.05
logFC_thrs_gene = 0
Uniprot_map_db = "/fsx/home/crivera/Data_Bases/UNIPROT_GeneInfo/hgnc_uniprot_mapping.txt"
workstr = "Proteomics"

Exps = Exps[Exps %in% c( "A-2025-0150-A", "A-2025-0150-B")]
for(exp in Exps){


    dir.create(file.path(OutputFolder, exp, "REPORTS"),  recursive = TRUE)
    List_contrasts <- read.csv(file.path(OutputFolder, exp,"CONTRAST_FILE_GEN", paste(exp, "constrasts.csv", sep = '_')))


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
                                        scinamicNum = exp,
                                        baseDir = file.path(OutputFolder, exp),
                                        DirDataForPipeline = DirDataForPipeline,
                                        DirPipeline = DirPipeline,
                                        Contrast_path = file.path(OutputFolder, exp,"CONTRAST_FILE_GEN", paste(exp, "constrasts.csv", sep = '_')),
                                        Chemo_path = Chemo_path,
                                        Experimental_design_path = "",
                                        Cell_Dict_path = Cell_Dict_path,
                                        Organism =  Organism,
                                        padj_thr_gene = padj_thr_gene,
                                        logFC_thrs_gene = logFC_thrs_gene,
                                        ReportNum = reportnum,
                                        Uniprot_map_db = Uniprot_map_db,
                                        Metadata_gr_path = file.path(baseDir_all, exp, paste(exp, "contrast_metadata.csv", sep = '_')),
                                        Contrasts_stats = file.path(baseDir_all, exp, paste(exp, "protein_abundance_contrast_stats.csv", sep = '_')),
                                        workstream = workstr
                          ),
                          clean = TRUE,
                          output_file = file.path(OutputFolder,exp, "REPORTS",paste(format(Sys.time(), '%y%m%d'), exp,"Report",reportnum, "DEG_ENRICH", '.html',sep = '_')))
    }
}



