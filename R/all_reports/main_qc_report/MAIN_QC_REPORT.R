library(rmarkdown)
library(tidyverse)
library(readxl)

source("/fsx/home/john/projects/empress/EXP000090/params_qc_groupeddonors.txt")


Name_folder_path =  file.path(baseDir, "REPORTS")
dir.create(Name_folder_path, recursive = TRUE)

# THe correlation plots are not being plotted using this function. We need to solve that using the childs or something else.

rmarkdown::render(input = paste(DirPipeline,"QC_REPORT.Rmd", sep = '/')  ,
                  params = list(cellline =  cellline,
                                scinamicNum = scinamicNum,
                                baseDir =  baseDir,
                                DirDataForPipeline = DirDataForPipeline,
                                DirPipeline = DirPipeline,
                                Metrics_path = Metrics_path,
                                Metadata_path = Metadata_path,
                                CountM_path =CountM_path,
                                SampleName_path = SampleName_path,
                                Experimental_design_path = Experimental_design_path,
                                Column_match = Column_match,
                                Organism_type = Organism_type,
                                QCNORM =  QCNORM,
                                TPM_path = TPM_path,
                                ControlName = ControlName,
                                Batch_variable = Batch_variable,
                                ExpDescription = Batch_variable,
                                Prev_perc =  Prev_perc,
                                MEMORY_MB = MEMORY_MB,
                                Corr_plot_features = Corr_plot_features,
                                CoarseCondition =  CoarseCondition,
                                Alias_path = Alias_path,
                                PRE_FILTER = PRE_FILTER),
                  clean = TRUE,
                  output_file = file.path(Name_folder_path,paste(format(Sys.time(), '%y-%m-%d'), cellline, "QC_REPORT", '.html',sep = '_')))


