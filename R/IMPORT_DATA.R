# 1. READS METADATA AND COUNT MATRIX
# 2. CHECKS COUNT MATRIX COLNAMES ARE CONTAINED WITHIN METADATA MATCH_FEATURE
# 3. ADDING METADATA FEATURES THAT WILL BE NEEDED LATER.

IMPORT_DATA = function(Metadata_path, Count_path, Match_feature, Output_file_path, SampleNamepath){

    library("openxlsx")
    library("readxl")



    dir.create(Output_file_path)
    setwd(Output_file_path)

    Metadata = read_xlsx(Metadata_path)
    Count = read.table(Input_matrix_path,  sep = '\t', header = TRUE, check.names=FALSE)

    # Check colnames of Count matrix are within the Metadata match_feature and visceversa

    if(any(!(colnames(Count) %in% Metadata[[Match_feature]])) == TRUE){

      print(!(colnames(Count) %in% Metadata[[Match_feature]]))
      stop("An element within the Count matrix is not in Metadata. Please check")

    }else if(any(!Metadata[[Match_feature]] %in% colnames(Count)) == TRUE){


     stop("An element within Metadata is not in Count matrix. Please check")
    }

    # Order samples correctly

    Count = Count[, match(Metadata[[Match_feature]], colnames(Count))]

    # Modify Treatment notation for convenience
    Metadata$Treatment_sci = Metadata$Treatment
    Metadata$Treatment = gsub("OL-00010", "", Metadata$Treatment)

    # Customizing sample names
    SampleNames <- read.csv(SampleNamepath, header = FALSE)
    SampleNames = SampleNames$V1
    Metadata$Sample_name = Metadata[[SampleNames[1]]]

    for(j in 2:length(SampleNames)){
        Metadata$Sample_name = paste(Metadata[["Sample_name"]],Metadata[[SampleNames[j]]], sep = "_")
    }

    # Replace NA by zero in Treatment concentration. Order concentrations
    if("Treatment_conc_uM" %in% colnames(Metadata)){

        Metadata$Treatment_conc_uM[Metadata$Treatment_conc_uM == "NA"] = 0
        Metadata$Treatment_conc_uM = factor(Metadata$Treatment_conc_uM, levels = sort(unique(as.numeric(Metadata$Treatment_conc_uM))))
    }else if("Treatment_concentration_mg_Kg" %in% colnames(Metadata)){

        Metadata$Treatment_concentration_mg_Kg[Metadata$Treatment_concentration_mg_Kg == "NA"] = 0
        Metadata$Treatment_concentration_mg_Kg = factor(Metadata$Treatment_concentration_mg_Kg, levels = sort(unique(as.numeric(Metadata$Treatment_concentration_mg_Kg))))

    }



    # Order collection time
    #Metadata$timeColl = factor(Metadata$`Rna_collection_time(hrs)`, levels = sort(as.numeric(unique(Metadata$`Rna_collection_time(hrs)`))))
    Metadata$timeColl = factor(Metadata$Rna_collection_time_hrs, levels = sort(as.numeric(unique(Metadata$Rna_collection_time_hrs))))

    # if(length(Suffix_del_name)!= 0){
    #
    #   print("Making sample name shorter")
    #   Metadata$Sample_name = gsub(Suffix_del_name, "", Metadata$Sample_id)
    #
    # }else{
    #
    #   Metadata$Sample_name = Metadata$Sample_id
    # }


    # Replacing column names of Count matrix to Sample_name

    colnames(Count) =  Metadata$Sample_name[match(colnames(Count), Metadata[[Match_feature]])]

    # Save matrix in output file

    dir.create("IMPORT_DATA")
    setwd("IMPORT_DATA")

    wb <- createWorkbook()
    addWorksheet(wb, "Count")
    addWorksheet(wb, "Metadata")

    writeData(wb, sheet = "Count", x = Count, rowNames = TRUE)
    writeData(wb, sheet = "Metadata", x = Metadata)

    saveWorkbook(wb, "Import_Data.xlsx", overwrite = TRUE)
    setwd(Output_file_path)

    return()


}
