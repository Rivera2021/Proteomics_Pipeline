
REPORT = function(Output_file_path, DirPipeline, DirDataForPipeline){

        library("rmdformats")
        library("rmarkdown")
        #library("readxl")
        #library("kableExtra")
        #library("dplyr")
        #library("ggplot2")
        #library("ggrepel")

        Dir_Output = Output_file_path
        DirDataForPipeline = DirDataForPipeline

        # Create Folder
        setwd(Output_file_path)
        Name_folder = "REPORT"
        dir.create(Name_folder)
        setwd(Name_folder)

        # Params
        Date = Sys.Date()

        Report_Name = paste(gsub('-', '_',Date), "Main_Report", sep = "_")
        rmarkdown::render(input = paste(DirPipeline, "Main_Report.Rmd", sep = '/'), params = list(Dir_Output = Dir_Output,DirDataForPipeline = DirDataForPipeline ),output_file = paste(Output_file_path,Name_folder,Report_Name, sep = "/"))

        # rmarkdown::render(input = paste(DirPipeline, "Main_Report.Rmd", sep = '/'),
        #                   params = list(Dir_Main = Dir_Main, Dir_Folders = Dir_Folders, Dir_Summary = Dir_Summary,
        #                                 List_File = List_subFold[[i]]), output_file = paste(Output_Folder,Report_Name, sep = "/"))

        setwd(Output_file_path)
}



