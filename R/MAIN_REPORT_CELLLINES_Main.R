MAIN_REPORT_CELLLINES_Main = function(Output_file_path, DirDataForPipeline,outliers_path, Cellline, DirPipeline, QCNORM, DEG_url_path, Check_Stim_Gene_path){


    library("rmdformats")
    library("rmarkdown")
    library("readxl")
    #library("kableExtra")
    library("dplyr")
    library("ggplot2")
    library("ggrepel")

    setwd(Output_file_path)
    # Create Folder
    Name_folder =  "MAIN_REPORTS"
    dir.create(Name_folder)
    setwd(Name_folder)



    rmarkdown::render(input = paste(DirPipeline,"MAIN_REPORT_CELLLINES.Rmd", sep = '/')  ,
                      params = list(Dir_Output = Output_file_path, DirDataForPipeline = DirDataForPipeline,  outliers_path =  outliers_path, Cellline = Cellline, QCNORM = QCNORM, DEG_url_path = DEG_url_path, Check_Stim_Gene_path = Check_Stim_Gene_path), output_file = paste(Output_file_path,Name_folder,paste("MAIN_REPORT", Cellline, '.html',sep = '_'),sep = "/"))

    setwd(Output_file_path)
}
