# QC- PRENORMALIZATION
# 1. ESTIMATED WHETHER THERE ARE OUTLIERS REGARDING NUMBER OF NON-ZERO READS PER SAMPLE OR NUMBER OF READS PER SAMPLE
# 2. PENDING: INCLUDE QC REGARDING CORRELATION WITHIN REPLICATES
# 3. PENDING: PLOTS SHOULD BE ADDED IN LOG SCALE


QC_PRENORMALIZATION_CELLLINE = function(Output_file_path,nMust = 2){
    library(tidyverse)
    library(dplyr)
    library("stringr")
    library("openxlsx")


    # Import data

    setwd(Output_file_path)
    Data_file = "./IMPORT_DATA/Import_Data.xlsx"
    Count = read.xlsx(xlsxFile = Data_file, sheet = "Count", rowNames= TRUE)
    Metadata = read.xlsx(xlsxFile = Data_file, sheet = "Metadata")

    # Create Folder

    Name_folder =  "QC_PRENORMALIZATION"
    dir.create(Name_folder)
    setwd(Name_folder)

    # Quality control on number of reads and  number of non-zero genes

    # Num Reads per sample
    Tot_Reads = colSums(Count)
    # Non-zero genes per sample
    Count_Genes = (Count !=0)
    Ngenes = colSums(Count_Genes)

    Metadata$ReadNum = Tot_Reads[match(Metadata$Sample_name,names(Tot_Reads))]
    Metadata$No_genes = Ngenes[match(Metadata$Sample_name,names(Ngenes))]

    # Convert NoReads and No_genes in log scale to estimate outliers in log scale

    Metadata$LogReadNum = log10(Metadata$ReadNum)
    Metadata$LogNo_genes = log10(Metadata$No_genes)

    # Thresholds for Number of reads and Number of non-zero genes
    IQR_reads = quantile(Metadata$LogReadNum,3/4) - quantile(Metadata$LogReadNum,1/4)
    Thrs_reads =  quantile(Metadata$LogReadNum,1/4) - nMust * IQR_reads

    IQR_genes = quantile(Metadata$LogNo_genes,3/4) - quantile(Metadata$LogNo_genes,1/4)
    Thrs_Nogenes =  quantile(Metadata$LogNo_genes,1/4) - nMust * IQR_genes

    Thrs_df = data.frame(Thrs_reads = Thrs_reads, Thrs_Nogenes = Thrs_Nogenes)

    # Find and save outliers

    MustDrop_List = c()
    Metadata_MustDrop = Metadata %>% filter(LogNo_genes < Thrs_Nogenes |  LogReadNum < Thrs_reads)
    MustDrop_list = Metadata_MustDrop$Sample_name


    Df = list()
    Df[["Dropped"]] = MustDrop_list
    write.xlsx(Df, file = 'Outliers.xlsx')

    # Plots

    pdf("QC_Prenormalization.pdf", width = 8, height = 6)

    # Boxplots of Logreads and Log non-zero genes with lower threshold
    p1 = ggplot(Metadata, aes(x = "",y = LogReadNum)) + geom_boxplot(outlier.shape = NA) +
        ylab("Log10 number of reads")+
        geom_jitter(position=position_jitter(w = 0.1, h = 0)) +
        geom_hline(yintercept = Thrs_reads, linetype = 'dashed', color = 'navy')+
        ggtitle(paste("Log10 of number of reads. Dashed: Q1 - ", as.character(nMust), "x IQR ",sep = " "))

    print(p1)

    p2 = ggplot(Metadata, aes(x = "",y = LogNo_genes)) + geom_boxplot(outlier.shape = NA) +
        ylab("Log10 number of non-zero genes")+
        geom_jitter(position=position_jitter(w = 0.1, h = 0)) +
        geom_hline(yintercept = Thrs_Nogenes, linetype = 'dashed', color = 'navy')+
        ggtitle(paste("Log10 number of non-zero genes. Dashed: Q1 - ", as.character(nMust), "x IQR ",sep = " "))

    print(p2)

    # Plot Number of reads Vs number of genes in log10 space with threshold for outliers



    p3 = ggplot(Metadata, aes(x=LogReadNum, y=LogNo_genes, group = Treatment)) +
        geom_point(aes(shape=timeColl, color=Treatment)) +
        xlab("Log10 of number of reads") +
        ylab("Log10 of number of non-zero genes") +
        xlim(min(Metadata$LogReadNum,Thrs_reads), max(Metadata$LogReadNum))+
        ylim(min(Metadata$LogNo_genes,Thrs_Nogenes), max(Metadata$LogNo_genes))+
        geom_hline(data = Thrs_df, aes(yintercept = Thrs_Nogenes), linetype="dashed")+
        geom_vline(data = Thrs_df, aes(xintercept = Thrs_reads), linetype="dashed")+
        ggtitle(paste("Number of reads vs Non-zero genes in Log10 space. Dashed: Q1 - ", as.character(nMust), "x IQR",sep = " "))

    print(p3)


    # Plot of Number reads per million vs Number of non-zero genes. In original space and show outliers selected
    Metadata$Out_QC_Prenorm = ""
    Metadata$Out_QC_Prenorm[Metadata$Sample_name %in% MustDrop_list] = Metadata$Sample_name[Metadata$Sample_name %in% MustDrop_list]

    p4 = ggplot(Metadata, aes(x=ReadNum, y=No_genes, group = Treatment)) +
        geom_point(aes(shape=timeColl, color=Treatment)) +
        geom_text(label=Metadata$Out_QC_Prenorm, nudge_x = 0.25, nudge_y = 0.25, check_overlap = F, size = 3) +
        xlab("Number of reads") +
        ylab("Number of non-zero genes")

    print(p4)

    dev.off()

    #save(MustDrop_list, OptDrop_List, Metadata_MustDrop,Metadata_OptDrop, file = "QC_Remove_lists.RData")

    # Save QC-ed data
    Metadata_QCPre = Metadata %>% filter(!Sample_name %in% MustDrop_list)
    Count_QCPre = Count[ ,match(Metadata_QCPre$Sample_name, colnames(Count))]

    Metadata = Metadata_QCPre
    Count = Count_QCPre
    wb <- createWorkbook()
    addWorksheet(wb, "Count")
    addWorksheet(wb, "Metadata")

    writeData(wb, sheet = "Count", x = Count, rowNames = TRUE)
    writeData(wb, sheet = "Metadata", x = Metadata)

    saveWorkbook(wb, "QC_Data.xlsx", overwrite = TRUE)

    setwd(Output_file_path)

    return()
}
