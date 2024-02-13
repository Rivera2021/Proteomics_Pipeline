
# QC- POSTNORMALIZATION
# 1. ESTIMATES OUTLIERS FOR THE ENTIRE SET OF SAMPLES
# 2. GETS RID OF SAMPLES WITH VERY LOW CORRELATION WITH ITS OWN REPLICATES

QC_POSTNORMALIZATION = function(Output_file_path, OUTLIER_FILTER = 'GENTLE_REP', RepFeatures_path){
  library(tidyverse)
  library(dplyr)
  library("stringr")
  library("PCAtools")
  library(BiocParallel)
  library(parallel)
  library("rstatix")
  library(umap)
  library(tidyverse)
  library("writexl")



  # Param for oulier
  nLow = 1
  nHigh = 3

  # Functions
  # Function to plot PCA with selected outliers
  Plot_UMAP_PCA = function(pc_proj_Sig,Metadata, OutlierToLabel , Folder_PCA){

        Name_folder =  Folder_PCA
        dir.create(Name_folder)
        setwd(Name_folder)


        # In order to uncomment this, we need to provide only significant components. Not worth to do it at the moment
        # Umap of outliers

        # umap_fit <- pc_proj_Sig %>% dplyr::select(where(is.numeric))  %>% scale() %>% umap()
        # umap_df <- umap_fit$layout %>% as.data.frame() %>% dplyr::rename(UMAP1="V1", UMAP2="V2")
        # umap_df$Sample_name = rownames(umap_df)
        # umap_df = umap_df %>% inner_join(Metadata, by="Sample_name")
        # umap_df$outlier = ""
        # umap_df$outlier[umap_df$Sample_name %in% OutlierToLabel] = umap_df$Sample_name[umap_df$Sample_name %in% OutlierToLabel]
        #
        #
        # # Plot UMAP
        #
        # P = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Treatment)) + geom_point()+labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot")+
        #     geom_text(label=umap_df$outlier, size = 3)
        #
        # pdf(file="UMAP_outliers.pdf", width = 12, height = 8)
        # print(P)
        # dev.off()

        # Plot PCA
        pc_proj_Sig_pca = pc_proj_Sig
        pc_proj_Sig_pca$Sample_name = rownames(pc_proj_Sig_pca)
        pc_proj_Sig_pca =  pc_proj_Sig_pca %>% inner_join(Metadata, by="Sample_name")
        pc_proj_Sig_pca$outlier = ""
        pc_proj_Sig_pca$outlier[pc_proj_Sig_pca$Sample_name %in% OutlierToLabel] = pc_proj_Sig_pca$Sample_name[pc_proj_Sig_pca$Sample_name %in% OutlierToLabel]

        P1 = ggplot(pc_proj_Sig_pca, aes(x = PC1, y = PC2, color = Treatment)) + geom_point()+labs(x = "PC1",y = "PC2",subtitle = "PCA plot") +
            geom_text(label=pc_proj_Sig_pca$outlier,
                      nudge_x = 0.25, nudge_y = 0.25,
                      check_overlap = F, size = 3)

        pdf(file="PC1_vs_PC2_outliers.pdf", width = 12, height = 8)
        print(P1)
        dev.off()

        P2 = ggplot(pc_proj_Sig_pca, aes(x = PC1, y = PC3, color = Treatment)) + geom_point()+labs(x = "PC1",y = "PC3",subtitle = "PCA plot") +
            geom_text(label=pc_proj_Sig_pca$outlier,
                      nudge_x = 0.25, nudge_y = 0.25,
                      check_overlap = F, size = 3)

        pdf(file="PC1_vs_PC3_outliers.pdf", width = 12, height = 8)
        print(P2)
        dev.off()

        P3 = ggplot(pc_proj_Sig_pca, aes(x = PC2, y = PC3, color = Treatment)) + geom_point()+labs(x = "PC2",y = "PC3",subtitle = "PCA plot") +
            geom_text(label=pc_proj_Sig_pca$outlier,
                      nudge_x = 0.25, nudge_y = 0.25,
                      check_overlap = F, size = 3)

        pdf(file="PC2_vs_PC3_outliers.pdf", width = 12, height = 8)
        print(P3)
        dev.off()

        P4 = ggplot(pc_proj_Sig_pca, aes(x = PC1, y = PC4, color = Treatment)) + geom_point()+labs(x = "PC1",y = "PC4",subtitle = "PCA plot") +
            geom_text(label=pc_proj_Sig_pca$outlier,
                      nudge_x = 0.25, nudge_y = 0.25,
                      check_overlap = F, size = 3)

        pdf(file="PC1_vs_PC4_outliers.pdf", width = 12, height = 8)
        print(P4)
        dev.off()

        setwd('..')
    }

  # Import data

  # Intro message
  print("Outlier detection using correlation within replicates")

  # Import data
  setwd(Output_file_path)
  #Data_file = load("~/InVivoTranscriptomics/Results/Results_231004/DESeq_Norm_NegBin_WOut/VST_ON_SVD_OFF_DESeq_res.RData")
  Data_file = "./DESEQ_NORM/DESeq_Norm.RData"
  load(Data_file)

  RepFeatures <- read.csv(RepFeatures_path, header = FALSE)


  # Create Folder

  Name_folder =  "QC_POSTNORMALIZATION"
  dir.create(Name_folder)
  setwd(Name_folder)

  # Create replicate feature. It is not necessarily the same as CoarseCondition used in normalization step

  RepFeatures = RepFeatures$V1
  Metadata$RepIdentify = Metadata[[RepFeatures[1]]]

  for(j in 2:length(RepFeatures)){
      Metadata$RepIdentify = paste(Metadata[["RepIdentify"]],Metadata[[RepFeatures[j]]], sep = "_")
  }


  # Data frames to save
  Corr_all = cor(NormCounts, method = "spearman")
  Corr_w_all = data.frame(Condition = character(), Corr_rep = numeric())
  for (cond in unique(Metadata$RepIdentify)){

      Metadata_cond =  Metadata %>% filter(RepIdentify == cond)
      Corr_w = Corr_all[Metadata_cond$Sample_name,Metadata_cond$Sample_name]
      Corr_wv = as.vector(Corr_w[upper.tri(Corr_w)])
      Df_temp = data.frame(Condition = cond, Corr_rep =  Corr_wv)
      Corr_w_all = rbind(Corr_w_all,Df_temp)


  }
  # Thresholds for Number of reads and Number of non-zero genes
  IQR_reads = quantile(Corr_w_all$Corr_rep,3/4) - quantile(Corr_w_all$Corr_rep,1/4)
  Thrs_Low = quantile(Corr_w_all$Corr_rep,1/4) - nLow * IQR_reads
  Thrs_High = quantile(Corr_w_all$Corr_rep,1/4) - nHigh * IQR_reads

  Thrs_reads =  data.frame(yintercept = c(Thrs_Low, Thrs_High), linetype = c("dashed","dotted"), name = c('Q1-IQR', 'Q1-3*IQR'))

  Plt = ggplot() +
      geom_point(data = Corr_w_all, aes(x = Condition, y = Corr_rep, colour = Condition))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
      geom_hline(data = Thrs_reads, aes(yintercept = yintercept, linetype = name))+
      scale_linetype_manual(
          "Thresholds", values = setNames(Thrs_reads$linetype, Thrs_reads$name)) +
      ggtitle(" Correlations within replicates per condition")

  pdf(file="Correlation_within_replicates.pdf", width = 12, height = 8)
  print(Plt)
  dev.off()

  if(OUTLIER_FILTER == 'GENTLE_REP'){

        # Spot entire bad conditions. Evaluate whether all correlations within same condition are lower than the low threshold

        Corr_w_all$Fail_ThrsLow = Corr_w_all$Corr_rep < Thrs_Low
        Corr_w_all_group = Corr_w_all %>% group_by(Condition) %>% summarize(Prod = prod(Fail_ThrsLow)) %>% filter(Prod == 1)

        if(nrow(Corr_w_all_group) != 0){
            Metadata_temp1 = Metadata %>% filter(RepIdentify == Corr_w_all_group$Condition)
            Bad_Samples_cond = Metadata_temp1$Sample_name
        } else{

            Bad_Samples_cond = c()
        }

        # Get rid of outlier samples. Evaluate whether the maximum correlation with other replicates is below the most exigent threshold
        Outlier_samples = c()
        for (sp in Metadata$Sample_name){

            Metadata_sp =  Metadata %>% filter(Sample_name == sp)
            # Select biological replicates
            Metadata_w = Metadata %>% filter(RepIdentify == Metadata_sp$RepIdentify)
            Corr_w = Corr_all[Metadata_sp$Sample_name,Metadata_w$Sample_name]
            Corr_w = Corr_w[-c(grep(Metadata_sp$Sample_name, names(Corr_w)))]
            if(max(Corr_w) < Thrs_High){
                Outlier_samples = c(Outlier_samples,sp)
            }

        }

        Outlier_selected = unique(c(Outlier_samples ,Bad_Samples_cond))

        # Metadata_filt = Metadata %>% filter(!Sample_name %in% Outlier_selected)
        # NormCounts_filt = NormCounts[,match(Metadata_filt$Sample_name, colnames(NormCounts))]



  }else if (OUTLIER_FILTER == 'NONE'){
      Bad_Samples_cond = c()
      Outlier_samples = c()
      Outlier_selected = c()
      # Metadata_filt = Metadata
      # NormCounts_filt = NormCounts

  }else{

      stop("Select a valid OUTLIER_FILTER")
  }

  # PCA plot all with outliers marked
  res.pca <- prcomp(t(NormCounts), scale. = TRUE)
  # The PC projections are stored in the "x" value of the prcomp object
  pc_proj <- as.data.frame(res.pca$x)

  Plot_UMAP_PCA(pc_proj_Sig = pc_proj,Metadata = Metadata,  OutlierToLabel = Outlier_selected , Folder_PCA = "PCA_WITH_OUTLIERS")

  # Save important data frames
  Outlier_DF = list()
  Outlier_DF[["Bad_Samples_cond"]] = data.frame(Bad_Samples_cond)
  Outlier_DF[["Outlier_samples"]] = data.frame(Outlier_samples)
  Outlier_DF[["Outlier_Selected"]] = data.frame(Outlier_selected)
  Outlier_DF[["Corr"]] = as.data.frame(Corr_all)

  writexl::write_xlsx(Outlier_DF, "Outlier_DF.xlsx")

  # Save list of outliers as .csv as well, so that it is accepted by DESEQ_NORM

  write.csv(data.frame(Outlier_selected), "Outliers_Selected.csv", row.names=FALSE)


  setwd(Output_file_path)


  return()
}















