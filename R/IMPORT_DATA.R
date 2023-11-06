# 1. READS METADATA AND COUNT MATRIX
# 2. CHECKS COUNT MATRIX COLNAMES ARE CONTAINED WITHIN METADATA MATCH_FEATURE
# 3. ADDING METADATA FEATURES THAT WILL BE NEEDED LATER.

IMPORT_DATA = function(Metadata_path, Count_path, Match_feature, Output_file_path, Suffix_del_name){
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

   Metadata = read_xlsx(Metadata_path)
   Count = Count[, match(Metadata[[Match_feature]], colnames(Count))]

   # Customizing sample names
   if(length(Suffix_del_name)!= 0){

     print("Making sample name shorter")
     Metadata$Sample_name = gsub(Suffix_del_name, "", Metadata$Sample_id)

   }else{

     Metadata$Sample_name = Metadata$Sample_id
   }


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
