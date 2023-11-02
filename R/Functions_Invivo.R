

Plot_PC_Invivo_V1 <- function(df, Metadata, Color_gg, Title){
  library(ggplot2)
  # library("factoextra")
  # df: Matrix nrow: # samples and ncol: # features
  # Metadata: Metadata with a column named Sample that coincides with the Sample names of df
  # sample.type.clrs: Uncomment scale_color_manual to manually add colors to these factor
  # Change groups to colored or shaped differently manually
  res.pca <- prcomp(df, scale. = TRUE)
  # The PC projections are stored in the "x" value of the prcomp object
  pc_proj <- res.pca$x
  
  ## Check these are really the projections
  #subDf =df
  #df_scaled = apply(subDf, 2, scale)
  #rownames(df_scaled) = rownames(subDf)
  #S = svd(df_scaled) 
  #Samp_Proj = df_scaled %*% S$v
  ##Yes! They are the same
  #Diff = max(Samp_Proj - pc_proj)
  
  # Lets continue building data frame to plot
  # convert to a tibble retaining the sample names as a new column
  pc_proj <- pc_proj %>% as_tibble(rownames = "Sample_name")
  # print the result
  pc_proj <- pc_proj %>% inner_join(Metadata, by="Sample_name")
  
  # Variance explained
  Percentage_Var = res.pca$sdev^2 *100/ sum(res.pca$sdev^2)
  
  pca_plot = list()
  Color_gg <- ensym(Color_gg)
  # PC1 vs PC2
  pca_plot[[1]] <- pc_proj %>% 
    # create the plot
    #ggplot(aes(x = PC1, y = PC2, color=!!Color_gg, label=id)) +
    ggplot(aes(x = PC1, y = PC2, color=!!Color_gg)) +
    geom_point(aes(shape=factor(IndicationOn)), alpha=0.9, size=2) +
    #geom_text_repel(size=3, max.overlaps = 60) +
    coord_fixed(ratio = 0.7) + 
    labs(title = Title, x = paste("PC1: (", format(Percentage_Var[1],digits=2),"%) ",sep = ""),
         y = paste("PC2: (", format(Percentage_Var[2],digits=2),"%) ",sep = ""))+
    theme_classic()
  # scale_shape_manual(values=c(21,24)) +
  #scale_color_manual(values = sample.type.clrs) 
  
  # PC1 vs PC3
  pca_plot[[2]] <- pc_proj %>% 
    # create the plot
    ggplot(aes(x = PC1, y = PC3, color=!!Color_gg)) +
    geom_point(aes(shape=factor(IndicationOn)), alpha=0.9, size=2) +
    #geom_text_repel(size=3, max.overlaps = 40) +
    coord_fixed(ratio = 0.7) + 
    labs(title = Title, x = paste("PC1: (", format(Percentage_Var[1],digits=2),"%) ",sep = ""),
         y = paste("PC3: (", format(Percentage_Var[3],digits=2),"%) ",sep = ""))+
    theme_classic()
  # scale_shape_manual(values=c(21,24)) +
  #scale_color_manual(values = cols) 
  
  # PC2 vs PC3
  pca_plot[[3]] <- pc_proj %>% 
    # create the plot
    ggplot(aes(x = PC2, y = PC3, color=!!Color_gg)) +
    geom_point(aes(shape=factor(IndicationOn)), alpha=0.9, size=2) +
    #geom_text_repel(size=3, max.overlaps = 20) +
    coord_fixed(ratio = 0.7) + 
    labs(title = Title, x = paste("PC2: (", format(Percentage_Var[2],digits=2),"%) ",sep = ""),
         y = paste("PC3: (", format(Percentage_Var[3],digits=2),"%) ",sep = ""))+
    theme_classic()
  # scale_shape_manual(values=c(21,24)) +
  #scale_color_manual(values = sample.type.clrs) 
  
  #pca_plot[[4]] <- fviz_eig(res.pca, ncp=10)
  
  return(pca_plot)
}

Plot_PC_Invivo_PerMol <- function(df, Metadata, Color_gg, Title){
  library(ggplot2)
  # library("factoextra")
  # df: Matrix nrow: # samples and ncol: # features
  # Metadata: Metadata with a column named Sample that coincides with the Sample names of df
  # sample.type.clrs: Uncomment scale_color_manual to manually add colors to these factor
  # Change groups to colored or shaped differently manually
  res.pca <- prcomp(df, scale. = TRUE)
  # The PC projections are stored in the "x" value of the prcomp object
  pc_proj <- res.pca$x
  
  ## Check these are really the projections
  #subDf =df
  #df_scaled = apply(subDf, 2, scale)
  #rownames(df_scaled) = rownames(subDf)
  #S = svd(df_scaled) 
  #Samp_Proj = df_scaled %*% S$v
  ##Yes! They are the same
  #Diff = max(Samp_Proj - pc_proj)
  
  # Lets continue building data frame to plot
  # convert to a tibble retaining the sample names as a new column
  pc_proj <- pc_proj %>% as_tibble(rownames = "Sample_name")
  # print the result
  pc_proj <- pc_proj %>% inner_join(Metadata, by="Sample_name")
  
  # Variance explained
  Percentage_Var = res.pca$sdev^2 *100/ sum(res.pca$sdev^2)
  
  pca_plot = list()
  Color_gg <- ensym(Color_gg)
  # PC1 vs PC2
  pca_plot[[1]] <- pc_proj %>% 
    # create the plot
    #ggplot(aes(x = PC1, y = PC2, color=!!Color_gg, label=id)) +
    ggplot(aes(x = PC1, y = PC2, color=!!Color_gg)) +
    geom_point(aes(shape=factor(timeColl)), alpha=0.9, size=2) +
    #geom_text_repel(size=3, max.overlaps = 60) +
    coord_fixed(ratio = 0.7) + 
    labs(title = Title, x = paste("PC1: (", format(Percentage_Var[1],digits=2),"%) ",sep = ""),
         y = paste("PC2: (", format(Percentage_Var[2],digits=2),"%) ",sep = ""))+
    theme_classic()
  # scale_shape_manual(values=c(21,24)) +
  #scale_color_manual(values = sample.type.clrs) 
  
  # PC1 vs PC3
  pca_plot[[2]] <- pc_proj %>% 
    # create the plot
    ggplot(aes(x = PC1, y = PC3, color=!!Color_gg)) +
    geom_point(aes(shape=factor(timeColl)), alpha=0.9, size=2) +
    #geom_text_repel(size=3, max.overlaps = 40) +
    coord_fixed(ratio = 0.7) + 
    labs(title = Title, x = paste("PC1: (", format(Percentage_Var[1],digits=2),"%) ",sep = ""),
         y = paste("PC3: (", format(Percentage_Var[3],digits=2),"%) ",sep = ""))+
    theme_classic()
  # scale_shape_manual(values=c(21,24)) +
  #scale_color_manual(values = cols) 
  
  # PC2 vs PC3
  pca_plot[[3]] <- pc_proj %>% 
    # create the plot
    ggplot(aes(x = PC2, y = PC3, color=!!Color_gg)) +
    geom_point(aes(shape=factor(timeColl)), alpha=0.9, size=2) +
    #geom_text_repel(size=3, max.overlaps = 20) +
    coord_fixed(ratio = 0.7) + 
    labs(title = Title, x = paste("PC2: (", format(Percentage_Var[2],digits=2),"%) ",sep = ""),
         y = paste("PC3: (", format(Percentage_Var[3],digits=2),"%) ",sep = ""))+
    theme_classic()
  # scale_shape_manual(values=c(21,24)) +
  #scale_color_manual(values = sample.type.clrs) 
  
  #pca_plot[[4]] <- fviz_eig(res.pca, ncp=10)
  
  return(pca_plot)
}

Volcano_plot_padj = function(res_Vol, padj_thr, logFC_thr, Title, maxOver,xLabel,yLabel){
  # Create Volcano plot
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
  # add a column of NAs
  res_Vol = res_Vol[!is.na(res_Vol$padj) ,]
  res_Vol$Diff_expr <- "NO"
  # if log2Foldchange > 5 and pvalue < 0.05, set as "UP" 
  res_Vol$Diff_expr[res_Vol$log2FoldChange >= logFC_thr & res_Vol$padj <= padj_thr] <- "UP"
  # if log2Foldchange < -5 and pvalue < 0.05, set as "DOWN"
  res_Vol$Diff_expr[res_Vol$log2FoldChange <= -logFC_thr & res_Vol$padj <= padj_thr] <- "DOWN"
  # Now write down the name of genes beside the points...
  # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  res_Vol$delabel <- NA
  res_Vol$delabel[res_Vol$Diff_expr != "NO"] <- rownames(res_Vol)[res_Vol$Diff_expr != "NO"]
  
  N_Sig = length(which(res_Vol$padj < 0.2)) 
  PV<-ggplot(data=res_Vol, aes(x=log2FoldChange, y=-log10(padj), col=Diff_expr, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel(size = 3) +
    scale_color_manual(values=c(DOWN = "steelblue", NO = "black", UP = "aquamarine3")) +
    geom_vline(xintercept=c(-logFC_thr,logFC_thr), col="navy", linetype="dotted") +
    geom_hline(yintercept=-log10(padj_thr), col="navy", linetype="dotted")+
    ggtitle(paste(Title, ", N_sig = ", format(N_Sig,digits=3), sep = ""))+
    xlab(xLabel)+
    ylab(yLabel)
  #setwd("~/Transcriptomic_Experiment/Data_from_nfcore/Volcano_plot_Exp1")
  
  
  return(PV)
  #setwd("~/Transcriptomic_Experiment")
}

Gene_Expression_plot = function(NormMatrix,Metadata, GeneList, ControlName){
  
  NormMatrix_filt = as.data.frame(NormMatrix[GeneList,])
  NormMatrix_filt$Gene = rownames(NormMatrix_filt)
  NormMatrix_Long = NormMatrix_filt %>% pivot_longer(!Gene, names_to = "Sample_name", values_to = "Norm_Expr") %>% as.data.frame()
  NormMatrix_Long$Gene = factor(NormMatrix_Long$Gene)
  NormMatrix_Long$time = Metadata$`Rna_collection_time(hrs)`[match(NormMatrix_Long$Sample_name, Metadata$Sample_name)]
  NormMatrix_Long$time = factor(NormMatrix_Long$time, levels = c(6,12,24,72))
  NormMatrix_Long$Treatment = Metadata$Treatment[match(NormMatrix_Long$Sample_name, Metadata$Sample_name)]
  NormMatrix_Long$Treatment = factor( NormMatrix_Long$Treatment, levels = c(ControlName, mol))
  NormMatrix_Long$Sample_name = NULL
  
  # Add very little noise to avoid T-test to complain
  NormMatrix_Long$Norm_Expr = NormMatrix_Long$Norm_Expr + rnorm(nrow(NormMatrix_Long),mean = 0,sd = 1e-7)
  
  # Generate pvalues.  step.increase = 0.06 seems to control of the pvalue test in the y axis.
  stat.test <- NormMatrix_Long %>% group_by(Gene, time) %>% t_test(Norm_Expr ~ Treatment, ref.group = "Vehicle")%>%add_significance()
  stat.test <- stat.test %>% add_xy_position(x = "time", dodge = 0.5, step.increase = 0.06) 
  #%>% rstatix::t_test(Norm_Expr ~ Treatment, ref.group = "Vehicle")
  # Plot 
  bxp <- ggboxplot(NormMatrix_Long, x = "time", y = "Norm_Expr", color = "Treatment", palette = "jco",
                   facet.by = "Gene", scales = "free", add = "dotplot")+
    stat_pvalue_manual(stat.test, label = "p.signif", size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.2)))
  
  
  return(bxp)
  
  
  
}



