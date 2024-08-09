# QC- PRENORMALIZATION
# 1. ESTIMATED WHETHER THERE ARE OUTLIERS REGARDING NUMBER OF NON-ZERO READS PER SAMPLE OR NUMBER OF READS PER SAMPLE
# 2. PENDING: INCLUDE QC REGARDING CORRELATION WITHIN REPLICATES
# 3. PENDING: PLOTS SHOULD BE ADDED IN LOG SCALE


QC_PRENORMALIZATION_CL = function(Output_file_path,nMust = 2,CellLine_Dict, CellLineDB_path, Gene_annotation, TPM_path){
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


    # Create Cell line dictionary for number of detected genes:

    # Make names in Metadata agree with the Data base
    CellLine_Dic <- read.csv(CellLine_Dict, header = TRUE)

    # EVENTUALLY THIS STEP WILL BE PROBABLY BE TRANSFER OUTSIDE OF THE PIPELINE AND WE WILL MANTAIN A DIRECTORY OUTSIDE
    # Getting a sense of expressed genes in cell lines
    #df_atlas <- readr::read_tsv(paste(params$DirDataForPipeline,"rna_celline_240330.tsv", sep = '/'))
    df_atlas <- readr::read_tsv(CellLineDB_path)
    colnames(df_atlas)<-str_replace_all(colnames(df_atlas), c(" " = "_"  ))

    Thr_cellline = data.frame(Cell_line = character(), Thrs_cellline = numeric())
    for(cellline in unique(Metadata$Cell_line)){

        celllinedb = CellLine_Dic$InDb[CellLine_Dic$Original == cellline]
        df_t = df_atlas %>% dplyr::filter(Cell_line == celllinedb & TPM > 0)
        if(Gene_annotation == 'Ensemble'){
          thr = ifelse(length(base::intersect(df_t$Gene, rownames(Count))) !=0, length(base::intersect(df_t$Gene, rownames(Count))),NA)
        }else if(Gene_annotation == 'Symbol'){
          thr = ifelse(length(base::intersect(df_t$Gene_name, rownames(Count))) !=0, length(base::intersect(df_t$Gene_name, rownames(Count))),NA)

        }
        Thr_cellline = Thr_cellline %>% add_row(Cell_line = cellline, Thrs_cellline = thr )
    }



    # Quality control on number of reads and  number of non-zero genes

    # Num Reads per sample
    Tot_Reads = colSums(Count)
    # Non-zero genes per sample
    Count_Genes = (Count !=0)
    Count_Genes = (Count >= 1)
    Ngenes = colSums(Count_Genes)

    Metadata$ReadNum = Tot_Reads[match(Metadata$Sample_name,names(Tot_Reads))]
    Metadata$detected_genes = Ngenes[match(Metadata$Sample_name,names(Ngenes))]

    # Convert NoReads and detected_genes in log scale to estimate outliers in log scale

    Metadata$LogReadNum = log10(Metadata$ReadNum)
    Metadata$Logdetected_genes = log10(Metadata$detected_genes)

    ## Thresholds for Number of reads and Number of non-zero genes PER CELL LINE

    IQR_reads = Metadata %>% group_by(Cell_line) %>% summarise(IQR = IQR(LogReadNum, na.rm = TRUE), first_quartile = quantile(LogReadNum, 0.25, na.rm = TRUE))
    IQR_reads$Thrs = IQR_reads$first_quartile - nMust * IQR_reads$IQR

    IQR_genes = Metadata %>% group_by(Cell_line) %>% summarise(IQR = IQR(Logdetected_genes, na.rm = TRUE), first_quartile = quantile(Logdetected_genes, 0.25, na.rm = TRUE))
    IQR_genes$Thrs = IQR_genes$first_quartile - nMust * IQR_genes$IQR





    # Find and save outliers

    MustDrop_List = c()
    Metadata$MustDrop = 'No'

    for(i in 1:nrow(Metadata)){

        cell_df_genes = IQR_genes %>% dplyr::filter(Cell_line == Metadata$Cell_line[i])
        cell_df_reads = IQR_reads %>% dplyr::filter(Cell_line == Metadata$Cell_line[i])

        if(Metadata$Logdetected_genes[i] < cell_df_genes$Thrs | Metadata$LogReadNum[i] < cell_df_reads$Thrs ){
            #print('here')
            Metadata$MustDrop[i] = 'Yes'
        }

    }


    MustDrop_list = Metadata %>% dplyr::filter(MustDrop == 'Yes')



    Df = MustDrop_list %>% dplyr::select(Sample_name, ID,Well,ReadNum, detected_genes,)
    write.xlsx(Df, file = paste(Output_file_path, "QC_PRENORMALIZATION",'Outliers.xlsx', sep = '/'))

    # Import TPM
    if(file.exists(TPM_path)){
        TPM_counts = read.table(TPM_path, header = TRUE)
        if('gene_id' %in% colnames(TPM_counts)){
            TPM_counts$gene_id = NULL
        }
        colnames(TPM_counts) = gsub('\\.','-',colnames(TPM_counts))
        # Generate Matrix with TPM expression per cell line

        TPM_counts_long = TPM_counts %>% pivot_longer(!gene_name, names_to = "Sample_name", values_to = "tpm")
        TPM_counts_long$Cell_line =  Metadata$Cell_line[match(TPM_counts_long$Sample_name, Metadata$Sample_name)]
        TPM_counts_expr = TPM_counts_long %>% dplyr::group_by(Cell_line,gene_name) %>% dplyr::summarize(median_tpm = median(tpm))

        # PLOTS TPMs
        TPM_counts_expr = TPM_counts_expr %>% dplyr::filter(median_tpm != 0 )
        saveRDS(TPM_counts_expr, file = file.path(Output_file_path, "QC_PRENORMALIZATION", "TPM_count_expr.RDS"))

    }else if(TPM_path == "3prime"){

        # 3' RNAseq we don't need to correct for gene length
        count_data2 = count_data %>% column_to_rownames(var="Gene")
        TPM_counts_expr = apply(count_data2, 2, function(x){x/sum(x)})



    }




    # PLOTS

    p1 = ggplot(Metadata, aes(x = "",y = LogReadNum)) + geom_boxplot(outlier.shape = NA) +
        facet_wrap(vars(Cell_line), ncol = min(length(unique(Metadata$Cell_line)),3))+
        ylab("Log10 number of reads")+
        geom_jitter(position=position_jitter(w = 0.1, h = 0)) +
        geom_hline(data = IQR_reads, aes(yintercept = Thrs),  linetype="dashed", color = 'navy')+
        ggtitle(paste("Log10 of number of reads. Dashed: Q1 - ", as.character(nMust), "x IQR ",sep = " ")) +
        theme(element_text(size=rel(14)))


    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_LogReads.pdf", sep = '/'), plot=p1, width = 8,height = 6, units = 'in')
    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_LogReads.png", sep = '/'), plot=p1, width = 8,height = 6, units = 'in')



    p2 = ggplot(Metadata, aes(x = "",y = Logdetected_genes)) + geom_boxplot(outlier.shape = NA) +
        facet_wrap(vars(Cell_line), ncol = min(length(unique(Metadata$Cell_line)),3))+
        ylab("Log10 number of non-zero genes")+
        geom_jitter(position=position_jitter(w = 0.1, h = 0)) +
        geom_hline(data = IQR_genes, aes(yintercept = Thrs),  linetype="dashed", color = 'navy')+
        ggtitle(paste("Log10 number of non-zero genes. Dashed: Q1 - ", as.character(nMust), "x IQR ",sep = " "))+
        theme(element_text(size=rel(14)))


    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_LogGenes.pdf", sep = '/'), plot=p2,width = 8,height = 6, units = 'in')
    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_LogGenes.png", sep = '/'), plot=p2,width = 8,height = 6, units = 'in')


    # Plot Number of reads Vs number of genes in log10 space with threshold for outliers


    p3 = ggplot(Metadata, aes(x=LogReadNum, y=Logdetected_genes, group = Treatment)) +
        geom_point(aes(shape=timeColl, color=Treatment)) +
        facet_wrap(vars(Cell_line), ncol = min(length(unique(Metadata$Cell_line)),3))+
        xlab("Log10 of number of reads") +
        ylab("Log10 of number of non-zero genes") +
        #xlim(min(Metadata$LogReadNum,Thrs_reads), max(Metadata$LogReadNum))+
        #ylim(min(Metadata$Logdetected_genes,Thrs_Nogenes), max(Metadata$Logdetected_genes))+
        geom_hline(data = IQR_genes, aes(yintercept = Thrs), linetype="dashed")+
        geom_hline(data = Thr_cellline, aes(yintercept = log10(Thrs_cellline)), linetype="solid")+
        geom_vline(data = IQR_reads, aes(xintercept = Thrs), linetype="dashed")+
        ggtitle(paste("Number of reads vs Non-zero genes in Log10 space. Dashed: Q1 - ", as.character(nMust), "x IQR, \n Solid: Number of detected genes from ATLAS",sep = " "))+
        theme(plot.title = element_text(size = 13))

    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_LogGenesVsLogReads.pdf", sep = '/'), plot=p3, width = 8,height = 6, units = 'in')
    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_LogGenesVsLogReads.png", sep = '/'), plot=p3, width = 8,height = 6, units = 'in')



    p5 = ggplot(Metadata, aes(x=LogReadNum, y=Logdetected_genes, group = Treatment)) +
        geom_point(aes(shape=timeColl, color=Replicate)) +
        facet_wrap(vars(Cell_line), ncol = min(length(unique(Metadata$Cell_line)),3))+
        xlab("Log10 of number of reads") +
        ylab("Log10 of number of non-zero genes") +
        #xlim(min(Metadata$LogReadNum,Thrs_reads), max(Metadata$LogReadNum))+
        #ylim(min(Metadata$Logdetected_genes,Thrs_Nogenes), max(Metadata$Logdetected_genes))+
        geom_hline(data = IQR_genes, aes(yintercept = Thrs), linetype="dashed")+
        geom_hline(data = Thr_cellline, aes(yintercept = log10(Thrs_cellline)), linetype="solid")+
        geom_vline(data = IQR_reads, aes(xintercept = Thrs), linetype="dashed")+
        ggtitle(paste("Number of reads vs Non-zero genes in Log10 space. Dashed: Q1 - ", as.character(nMust), "x IQR, \n Solid: Number of detected genes from ATLAS",sep = " "))+
        theme(plot.title = element_text(size = 13))
    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_LogGenesVsLogReads_V2.pdf", sep = '/'), plot=p5, width = 8,height = 6, units = 'in')
    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_LogGenesVsLogReads_V2.png", sep = '/'), plot=p5, width = 8,height = 6, units = 'in')



    p6 = ggplot(Metadata, aes(x=LogReadNum, y=Logdetected_genes, group = Treatment)) +
        geom_point(aes(shape=timeColl, color=RQN)) +
        facet_wrap(vars(Cell_line), ncol = min(length(unique(Metadata$Cell_line)),3))+
        xlab("Log10 of number of reads") +
        ylab("Log10 of number of non-zero genes") +
        #xlim(min(Metadata$LogReadNum,Thrs_reads), max(Metadata$LogReadNum))+
        #ylim(min(Metadata$Logdetected_genes,Thrs_Nogenes), max(Metadata$Logdetected_genes))+
        geom_hline(data = IQR_genes, aes(yintercept = Thrs), linetype="dashed")+
        geom_hline(data = Thr_cellline, aes(yintercept = log10(Thrs_cellline)), linetype="solid")+
        geom_vline(data = IQR_reads, aes(xintercept = Thrs), linetype="dashed")+
        ggtitle(paste("Number of reads vs Non-zero genes in Log10 space. Dashed: Q1 - ", as.character(nMust), "x IQR, \n Solid: Number of detected genes from ATLAS",sep = " "))+
        theme(plot.title = element_text(size = 13))

    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_LogGenesVsLogReads_V3.pdf", sep = '/'), plot=p6, width = 8,height = 6, units = 'in')
    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_LogGenesVsLogReads_V3.png", sep = '/'), plot=p6, width = 8,height = 6, units = 'in')




    # Plot of Number reads per million vs Number of non-zero genes. In original space and show outliers selected
    Metadata$Out_QC_Prenorm = ""
    Metadata$Out_QC_Prenorm[Metadata$Sample_name %in% Df$Sample_name] = Metadata$Sample_name[Metadata$Sample_name %in% Df$Sample_name]



    p4 = ggplot(Metadata, aes(x=ReadNum, y=detected_genes, group = Treatment)) +
        geom_point(aes(shape=timeColl, color=Treatment)) +
        facet_wrap(vars(Cell_line), ncol = min(length(unique(Metadata$Cell_line)),3), scales = 'free') +
        geom_hline(data = IQR_genes, aes(yintercept = 10^(Thrs)), linetype="dashed")+
        geom_hline(data = Thr_cellline, aes(yintercept = Thrs_cellline), linetype="solid")+
        geom_vline(data = IQR_reads, aes(xintercept = 10^(Thrs)), linetype="dashed")+
        #geom_text(label=Metadata$Out_QC_Prenorm, nudge_x = 0, nudge_y = 0, check_overlap = T, size = 2) +
        xlab("Number of reads") +
        ylab("Number of non-zero genes")+
        ggtitle(paste("Number of reads vs Non-zero genese. Dashed: Q1 - ", as.character(nMust), "x IQR,  \n Solid: Number of detected genes from ATLAS",sep = " "))+
        theme(plot.title = element_text(size = 15))

    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_GenesVsReads.pdf", sep = '/'), plot=p4, width = 8,height = 6, units = 'in')
    ggsave(filename=paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Prenorm_GenesVsReads.png", sep = '/'), plot=p4, width = 8,height = 6, units = 'in')


    # Save QC-ed data
    # Metadata_QCPre = Metadata %>% filter(!Sample_name %in% Df$Sample_name)
    # Count_QCPre = Count[ ,match(Metadata_QCPre$Sample_name, colnames(Count))]
    #
    # Metadata = Metadata_QCPre
    # Count = Count_QCPre
    wb <- createWorkbook()
    addWorksheet(wb, "Count")
    addWorksheet(wb, "Metadata")

    writeData(wb, sheet = "Count", x = Count, rowNames = TRUE)
    writeData(wb, sheet = "Metadata", x = Metadata)

    saveWorkbook(wb, paste(Output_file_path, "QC_PRENORMALIZATION", "QC_Data.xlsx", sep = '/'), overwrite = TRUE)


    setwd(Output_file_path)

    return()
}

