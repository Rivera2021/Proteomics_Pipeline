
# QC- PRENORMALIZATION
# 1. ESTIMATED WHETHER THERE ARE OUTLIERS REGARDING NUMBER OF NON-ZERO READS PER SAMPLE OR NUMBER OF READS PER SAMPLE
# 2. PENDING: INCLUDE QC REGARDING CORRELATION WITHIN REPLICATES
# 3. PENDING: PLOTS SHOULD BE ADDED IN LOG SCALE


QC_PRENORMALIZATION = function(Output_file_path, nOpt, nMust){
    library(tidyverse)
    library(dplyr)
    library("stringr")


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


    n = nOpt
    nmax = nMust

    # Thresholds for Number of reads and Number of non-zero genes
    Thrs_Reads = data.frame(Mean = mean(Metadata$ReadNum), Sd = sd(Metadata$ReadNum))
    Thrs_Reads$Thrs = (Thrs_Reads$Mean - n*Thrs_Reads$Sd)/1e6
    Thrs_Reads$ThrsMust = (Thrs_Reads$Mean - nmax*Thrs_Reads$Sd)/1e6

    Thrs_NoGenes = data.frame(Mean = mean(Metadata$No_genes), Sd = sd(Metadata$No_genes))
    Thrs_NoGenes$Thrs = (Thrs_NoGenes$Mean - n*Thrs_NoGenes$Sd)
    Thrs_NoGenes$ThrsMust = (Thrs_NoGenes$Mean - nmax*Thrs_NoGenes$Sd)


    pdf("QC_Prenormalization.pdf", width = 8, height = 6)

    p1 = hist(Metadata$ReadNum, 20, main = "Hist reads/sample and 2 SD", xlab = "No Reads")
    abline(v = Thrs_Reads$Mean + n*Thrs_Reads$Sd, col="red")
    abline(v = Thrs_Reads$Mean - n*Thrs_Reads$Sd, col="red")
    print(p1)

    p2 = hist(Metadata$No_genes, 20, main = "Hist non-zero genes/sample and 2 SD", xlab = "Non-zero genes")
    abline(v = Thrs_NoGenes$Mean + n*Thrs_NoGenes$Sd, col="red")
    abline(v = Thrs_NoGenes$Mean - n*Thrs_NoGenes$Sd, col="red")
    print(p2)



    p4 = ggplot(Metadata, aes(x=ReadNum/1e6, y=No_genes, color=as.factor(Treatment))) + geom_point() +
        xlab("Reads Per Million") +
        ylab("Non-zero genes") +
        xlim(min(Metadata$ReadNum/1e6,Thrs_Reads$ThrsMust), max(Metadata$ReadNum/1e6))+
        ylim(min(Metadata$No_genes,Thrs_NoGenes$ThrsMust), max(Metadata$No_genes))+
        geom_hline(data = Thrs_NoGenes, aes(yintercept = ThrsMust), linetype="dashed")+
        geom_vline(data = Thrs_Reads, aes(xintercept = ThrsMust), linetype="dashed")+
        ggtitle(paste("RPM vs Non-zero genes. Dashed:", as.character(nmax), "SD",sep = " "))


    print(p4)



    Metadata$timeColl = factor(Metadata$`Rna_collection_time(hrs)`, levels = sort(as.numeric(unique(Metadata$`Rna_collection_time(hrs)`))))
    p5 = ggplot(Metadata, aes(x=ReadNum/1e6, y=No_genes, group = Treatment)) + geom_point(aes(shape=timeColl, color=Treatment)) +
        xlab("Reads Per Million") +
        ylab("Non-zero genes") +
        xlim(min(Metadata$ReadNum/1e6,Thrs_Reads$ThrsMust), max(Metadata$ReadNum/1e6))+
        ylim(min(Metadata$No_genes,Thrs_NoGenes$ThrsMust), max(Metadata$No_genes))+
        geom_hline(data = Thrs_NoGenes, aes(yintercept = ThrsMust), linetype="dashed")+
        geom_vline(data = Thrs_Reads, aes(xintercept = ThrsMust), linetype="dashed")+
        ggtitle(paste("RPM vs Non-zero genes. Dashed:", as.character(nmax), "SD",sep = " "))
    print(p5)

    # Rna concentration has some NA
    Metadata$Rna_concentration = Metadata$Rna_concentration
    Metadata_fil = Metadata %>% filter(!is.na(Rna_concentration))
    Metadata_fil$Rna_concentration = as.numeric(Metadata_fil$Rna_concentration)

    p6 = ggplot(Metadata_fil, aes(x=timeColl, y=Rna_concentration, color = timeColl)) + geom_boxplot() +
        xlab("time collection") +
        ylab("Rna_concentration") +
        ggtitle(paste("Time vs RNA concentration",sep = " "))+
        geom_jitter(shape=16, position=position_jitter(0.2))
    print(p6)



    dev.off()

    # Get a list of must drop and a list of possible candidates to drop out.
    MustDrop_List = c()
    OptDrop_List = c()
    Metadata_Exps_cell_Must = list()
    Metadata_Exps_cell_opt = list()
    Metadata_MustDrop = Metadata %>% filter(No_genes < Thrs_NoGenes$ThrsMust |  ReadNum < Thrs_Reads$ThrsMust * 1e6)
    Metadata_OptDrop = Metadata %>% filter(No_genes < Thrs_NoGenes$Thrs |  ReadNum < Thrs_Reads$Thrs * 1e6)
    MustDrop_list = Metadata_MustDrop$Sample_name
    OptDrop_List = Metadata_OptDrop$Sample_name
    Df = list()
    Df[["Dropped"]] = MustDrop_list
    Df[["OptDrop"]] = OptDrop_List

    write.xlsx(Df, file = 'Outliers.xlsx')
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

