# PRE-FILTERING
# 1. REDUCE NUMBER OF GENES EITHER BY PREVALENCE OR BY MIXTURE OF NEGATIVE BINOMIAL

PRE_FILTERING = function(Output_file_path, Prev_perc, PRE_FILTER, CoarseConditions){

    library("openxlsx")
    library("dplyr")

    # Functions
    Filtered_Prevalence = function(Matrix, Prev_perc ){
        # nrow of Matrix: No of features
        # ncol of Matrix: No of samples

        Matrix_prev = as.matrix(Matrix %>% mutate_if(is.numeric, ~1 * (. != 0)))
        Indx = which(rowSums(Matrix_prev) >= ncol(Matrix)*Prev_perc)
        Matrix_Filt = Matrix[Indx,]
        return(Matrix_Filt)

    }

    # Import data

    setwd(Output_file_path)

    Data_file = "./QC_PRENORMALIZATION/QC_Data.xlsx"
    Count = read.xlsx(xlsxFile = Data_file, sheet = "Count", rowNames= TRUE)
    Metadata = read.xlsx(xlsxFile = Data_file, sheet = "Metadata")

    # # Get rid of outliers
    #
    # Outliers = read.xlsx(xlsxFile = Outliers_QC_prenorm_path)
    # Metadata = Metadata %>% filter(!Sample_name %in% Outliers$Sample_name)
    # Count = Count[, match(Metadata$Sample_name, colnames(Count))]

    # Convert to integer
    Count = Count %>% mutate_if(is.numeric, round)

    ##Create coarse Condition feature and useful Metadata
    Metadata$CoarseCondition = Metadata[[CoarseConditions[1]]]

    for(j in 2:length(CoarseConditions)){
        Metadata$CoarseCondition = paste(Metadata[["CoarseCondition"]],Metadata[[CoarseConditions[j]]], sep = "_")
    }


    # Pre-filter
    if(PRE_FILTER =='PREV'){

       #print("Pre-filtering using prevalence")
       Count_Filt = Filtered_Prevalence(Count, Prev_perc )
       Metadata = Metadata[match(colnames(Count_Filt), Metadata$Sample_name),]
    }else if(PRE_FILTER == 'LOW_EXPR'){
        Cond_Num = Metadata %>% group_by(CoarseCondition) %>% summarize(group_size = n())
        smallestGroupSize <- min(Cond_Num$group_size)
        keep <- rowSums(Count >= 10) >= smallestGroupSize
        Count_Filt <- Count[keep,]
        Metadata = Metadata[match(colnames(Count_Filt), Metadata$Sample_name),]
    }
    else if(PRE_FILTER =='MNB'){

       #print("Pre-filtering using mixture of negative binomial")





    }



   # To get the set of expressed genes use for chemoproteomics we lower the threshold a bit
    Cond_Num = Metadata %>% group_by(CoarseCondition) %>% summarize(group_size = n())
    smallestGroupSize <- min(Cond_Num$group_size)
    keep_l <- rowSums(Count >= 1) >= smallestGroupSize
    Count_expr <- Count[keep_l,]


   # Create Folder

   Name_folder =  "PRE_FILTERING"
   dir.create(Name_folder)
   setwd(Name_folder)

   # save matrix and metadata
   Count = Count_Filt
   wb <- createWorkbook()
   addWorksheet(wb, "Count")
   addWorksheet(wb, "Metadata")
   writeData(wb, sheet = "Count", x = Count, rowNames = TRUE)
   writeData(wb, sheet = "Metadata", x = Metadata)
   saveWorkbook(wb, "Prefilter_Data.xlsx", overwrite = TRUE)

   # save count metrix for expressed genes
   wb <- createWorkbook()
   addWorksheet(wb, "Count_expr")
   writeData(wb, sheet = "Count_expr", x = Count_expr, rowNames = TRUE)
   saveWorkbook(wb, "Expressed_genes.xlsx", overwrite = TRUE)

   setwd(Output_file_path)

   return()

}



