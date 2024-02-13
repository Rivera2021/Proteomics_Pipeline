
REPORT = function(Output_file_path, DirPipeline ){

        library("rmdformats")
        library("rmarkdown")
        #library("readxl")
        #library("kableExtra")
        #library("dplyr")
        #library("ggplot2")
        #library("ggrepel")

        # Dir = "/fsx/home/crivera/Novogene/usftp21.novogene.com_Exp5"
        # Dir_Folders = "/fsx/home/crivera/Novogene/usftp21.novogene.com_Exp5/Results_HBEC5i/DeSeq_Resu_WOut_NfCore_Pvalue"
        # Dir_Summary =  "/fsx/home/crivera/Novogene/usftp21.novogene.com_Exp5/Results_HBEC5i/Summary_Plots_Nfcore"
        # Cell_type = "HBEC5i_LPS"


        # Create Folder
        setwd(Output_file_path)
        Name_folder =  "REPORT"
        dir.create(Name_folder)
        setwd(Name_folder)

        # Params
        Date = Sys.Date()

        Report_Name = paste(gsub('-', '_',Date), "Main_Report", sep = "_")
        rmarkdown::render(input = paste(DirPipeline, "Main_Report.Rmd", sep = '/'), output_file = paste(Output_file_path,Name_folder,Report_Name, sep = "/"))

        # rmarkdown::render(input = paste(DirPipeline, "Main_Report.Rmd", sep = '/'),
        #                   params = list(Dir_Main = Dir_Main, Dir_Folders = Dir_Folders, Dir_Summary = Dir_Summary,
        #                                 List_File = List_subFold[[i]]), output_file = paste(Output_Folder,Report_Name, sep = "/"))



}



