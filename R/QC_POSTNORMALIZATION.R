
# QC- POSTNORMALIZATION
# 1. ESTIMATES OUTLIERS FOR THE ENTIRE SET OF SAMPLES
# 2. GETS RID OF SAMPLES WITH VERY LOW CORRELATION WITH ITS OWN REPLICATES



QC_POSTNORMALIZATION = function(Output_file_path, OUTLIER_FILTER = 'GENTLE', RepFeatures_path){
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

  # OUTLIER_FILTER : GENTLE (Not all outliers predicted from Mahalanobis distance are removed. Only those that ALSO have a correlation
  # within replicate less than 0.8. After removing these outliers, correlation within replicate is recomputed, and those samples with average
  # correlation within replicate less than 0.8 are removed. ) STRONG (All candidate outliers from Mahalanobis distance are removed. After
  # removing these samples correlation within replicate is compute. All samples with average correlation less than 0.8 get removed).NONE
  # (No outlier detection is done)
  # The removed samples should be the union of Outlier_over and Outlier_rep within the Outlier_DF


  # Functions
    # Function to find Mahalanobis distance from of a point from the set of points, but controls for possible
    # outliers within the set of points selected to estimate the multi-gaussian distribution.
    mahalanobis.trim_v2 = function(y,x) {
        n = ncol(x)
        outm = apply(x,2,function(obj) is_outlier(obj,coef=3.0) )
        nout = apply(outm,1,sum)
        ind = nout==0
        xm = x[ind,]
        cv = matrix(0,n,n)
        mu = apply(xm,2,function(obj) mean(obj))
        for (i in 1:(n-1)) {
            for (j in (i+1):n) {
                cv[i,j] = mean( (xm[,i]-mu[i])*(xm[,j]-mu[j]) )
                cv[j,i] = cv[i,j]
            }
        }
        ### trimmed variances
        for (i in 1:n) {cv[i,i] =  mean( (xm[,i]-mu[i])^2 ) }
        ### assemble
        ans =  t(y-mu) %*% solve(cv) %*% (y-mu)
        ans = sqrt(ans)
        return(ans)
    }
    # Function to plot PCA with selected outliers
    Plot_UMAP_PCA = function(pc_proj_Sig,Metadata, OutlierToLabel , Folder_PCA){

        Name_folder =  Folder_PCA
        dir.create(Name_folder)
        setwd(Name_folder)

        # Umap of outliers
        umap_fit <- pc_proj_Sig %>% dplyr::select(where(is.numeric))  %>% scale() %>% umap()
        umap_df <- umap_fit$layout %>% as.data.frame() %>% dplyr::rename(UMAP1="V1", UMAP2="V2")
        umap_df$Sample_name = rownames(umap_df)
        umap_df = umap_df %>% inner_join(Metadata, by="Sample_name")
        umap_df$outlier = ""
        umap_df$outlier[umap_df$Sample_name %in% OutlierToLabel] = umap_df$Sample_name[umap_df$Sample_name %in% OutlierToLabel]


        # Plot UMAP

        P = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Treatment)) + geom_point()+labs(x = "UMAP1",y = "UMAP2",subtitle = "UMAP plot")+
            geom_text(label=umap_df$outlier, size = 3)

        pdf(file="UMAP_outliers.pdf", width = 12, height = 8)
        print(P)
        dev.off()

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
  print("Outlier detection using all samples...")

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

  # Data frames to save
  Outlier_DF = list()
  MH_pval = 0.001
  MH_Conf = 1-MH_pval
  Thrs_Corr = 0.8

  # MAHALANOBIS--------------------------------------------------------------------
  # Mahalanobis distance is estimated in PCA space.

  # Get significant PCA
  hornRes = parallelPCA(
      NormCounts,
      max.rank = min(nrow(NormCounts), ncol(NormCounts)),
      niters = 200,
      center = TRUE,
      scale=TRUE,
      threshold = 0.1,
      transposed = FALSE,
      BPPARAM = DoparParam()
  )

  print(paste("Number of significant PCs is ",as.character(hornRes$n), sep = ' ' ))

  pc_proj = hornRes$original$rotated
  pc_proj_Sig = pc_proj[,1:hornRes$n]

  # Mahalanobis distance method developed by Tony to estimate outliers. This method lets us figure out our critical threshold. Instead of
  # using Chi-square directly we estimate the critical value from simulated data. Values higher than this have pvalue <0.05
  ###Ensemble of 1000. Using the dimensions of out problem
  xs = rep(NA,1000)
  for (i in 1:1000) {
      x = matrix(rnorm(nrow(pc_proj_Sig)*ncol(pc_proj_Sig)), nrow=nrow(pc_proj_Sig), ncol=ncol(pc_proj_Sig))
      #x[1,] = c(10,0,0)
      y = rnorm(ncol(pc_proj_Sig))
      xs[i] = mahalanobis.trim_v2(y,x)^2
      #xs[i] = mahalanobis.trim_v1(y,x,0.1)^2
  }

  MH_crt = quantile(xs, MH_Conf)
  #qchisq(MH_Conf, ncol(pc_proj_Sig))

  # Estimate Mahalanobis distances for our data set

  MH_dist = c()
  for(j in 1:nrow(pc_proj_Sig)){

      y = as.numeric(pc_proj_Sig[j,])
      x = as.matrix(pc_proj_Sig[-j,])

      MH_dist[rownames(pc_proj_Sig)[j]] = mahalanobis.trim_v2(y,x)^2
      #MH_dist[rownames(pc_proj_Sig)[j]] = mahalanobis.trim_v1(y,x, 0.1)^2

  }

  Outlier_MH = MH_dist[MH_dist>MH_crt]
  Outlier_MH = Outlier_MH[order(Outlier_MH, decreasing = TRUE)]

  MH_dist_df = data.frame("Sample_Name" = names(MH_dist), "MH_dist" = MH_dist)


  # PLOTS MAHALANOBIS--------------------------------------------------------------------------
  # Plot Mahalanobis distance and threshold

  pdf(file="Mahalanobis_screeplot.pdf", width = 12, height = 8)
  plot(sort(MH_dist, decreasing = TRUE), xlab = "Rank", ylab = 'Mahalanobis distance')
  abline(h = MH_crt,  lty = 2)
  legend("topright", legend=c("MH_dist", "MH_threshold"),
          lty=1:2, cex=0.8,
         title="Line types")
  dev.off()

  Plot_UMAP_PCA(pc_proj_Sig,Metadata, OutlierToLabel = names(Outlier_MH) , Folder_PCA = "MAHALANOBIS_OUTLIERS_PCA")

   # CORRELATION WITHIN AND BETWEEN CONDITIONS---------------------------------------------------------

   print("Removing samples with low correlation within replicates ...")
   Corr_all = cor(NormCounts)
   Outlier_WRep = c()

   # Create replicate feature. It is not necessarily the same as CoarseCondition used in normalization step

   RepFeatures = RepFeatures$V1
   Metadata$RepIdentify = Metadata[[RepFeatures[1]]]

   for(j in 2:length(RepFeatures)){
       Metadata$RepIdentify = paste(Metadata[["RepIdentify"]],Metadata[[RepFeatures[j]]], sep = "_")
   }

   for (sp in Metadata$Sample_name){

       Metadata_sp =  Metadata %>% filter(Sample_name == sp)
       # Select biological replicates
       Metadata_w = Metadata %>% filter(RepIdentify == Metadata_sp$RepIdentify)
       Corr_w = Corr_all[Metadata_sp$Sample_name,Metadata_w$Sample_name]
       Corr_w = Corr_w[-c(grep(Metadata_sp$Sample_name, names(Corr_w)))]
       if (any(Corr_w < Thrs_Corr)){
           Outlier_WRep = c(Outlier_WRep,Corr_w[Corr_w < Thrs_Corr])
       }

   }


   if(OUTLIER_FILTER == 'GENTLE'){

       Outlier_Over = intersect(names(Outlier_MH), names(Outlier_WRep))
       # There is a set of points that have to be eliminated when the MH_dist show extreme values, even if the correlation within replicate points
       # is not less than 0.8
       Outlier_outlier_MH = intersect(names(MH_dist)[is_outlier(MH_dist,coef=3.0)],  names(Outlier_MH))
       Outlier_Over = unique(c(Outlier_Over,Outlier_outlier_MH))
       # Remove Outlier_over
       Metadata_filt = Metadata %>% filter(!Sample_name %in% Outlier_Over)
       NormCounts_filt = NormCounts[,match(Metadata_filt$Sample_name, colnames(NormCounts))]

       # Recompute correlations within replicate
       Corr_all_filt = cor(NormCounts_filt)
       Outlier_WRep_V2 = c()
       Mean_within = c()
       for (sp in Metadata_filt$Sample_name){

           Metadata_sp =  Metadata_filt %>% filter(Sample_name == sp)
           # Select biological replicates
           Metadata_w = Metadata_filt %>% filter(RepIdentify == Metadata_sp$RepIdentify)
           Corr_w = Corr_all[Metadata_sp$Sample_name,Metadata_w$Sample_name]
           Corr_w = Corr_w[-c(grep(Metadata_sp$Sample_name, names(Corr_w)))]
           if(length(Corr_w) == 0 ){
               next
           }
           Mean_within[sp] = mean(Corr_w)

           if (any(Mean_within[sp] < Thrs_Corr)){

               Outlier_WRep_V2 = c(Outlier_WRep_V2,Mean_within[sp])
           }

       }

       # Remove those last samples that have low correlation within replicate

       Metadata_filt = Metadata_filt %>% filter(!Sample_name %in% names(Outlier_WRep_V2))
       NormCounts_filt = NormCounts_filt[,match(Metadata_filt$Sample_name, colnames(NormCounts_filt))]



   }else if (OUTLIER_FILTER == 'STRONG'){

       Outlier_Over = names(Outlier_MH)
       # Remove Outlier_over
       Metadata_filt = Metadata %>% filter(!Sample_name %in% Outlier_Over)
       NormCounts_filt = NormCounts[,match(Metadata_filt$Sample_name, colnames(NormCounts))]


       # Compute within replicate correlation
       Corr_all_filt = cor(NormCounts_filt)
       Outlier_WRep_V2 = c()
       Mean_within = c()
       for (sp in Metadata_filt$Sample_name){

           Metadata_sp =  Metadata_filt %>% filter(Sample_name == sp)
           # Select biological replicates
           Metadata_w = Metadata_filt %>% filter(RepIdentify == Metadata_sp$RepIdentify)
           Corr_w = Corr_all[Metadata_sp$Sample_name,Metadata_w$Sample_name]
           Corr_w = Corr_w[-c(grep(Metadata_sp$Sample_name, names(Corr_w)))]
           if(length(Corr_w) == 0 ){
               next
           }
           Mean_within[sp] = mean(Corr_w)

           if (any(Mean_within[sp] < Thrs_Corr)){

               Outlier_WRep_V2 = c(Outlier_WRep_V2,Mean_within[sp])
           }
       }
           # Remove those last samples that have low correlation within replicate

           Metadata_filt = Metadata_filt %>% filter(!Sample_name %in% names(Outlier_WRep_V2))
           NormCounts_filt = NormCounts_filt[,match(Metadata_filt$Sample_name, colnames(NormCounts_filt))]

    }else if (OUTLIER_FILTER == 'NONE'){
           Outlier_Over = c()
           Outlier_WRep_V2 = c()
           Metadata_filt = Metadata
           NormCounts_filt = NormCounts

    }else{

        stop("Select a valid OUTLIER_FILTER")
    }

   Outlier_selected = c(Outlier_Over, names(Outlier_WRep_V2))
   Plot_UMAP_PCA(pc_proj_Sig,Metadata, OutlierToLabel = Outlier_selected, Folder_PCA = "SELECTED_OUTLIERS_PCA")

   # Save important data frames

   Outlier_DF[["MH_dist_df"]] = MH_dist_df
   Outlier_DF[["MH_crt"]] = data.frame(MH_crt)
   Outlier_DF[["Corr"]] = as.data.frame(Corr_all)
   Outlier_DF[["Outlier_Over"]] = data.frame(Outlier_Over)
   Outlier_DF[["Outlier_RepV2"]] = data.frame(Outlier_WRep_V2)
   Outlier_DF[["Outlier_Selected"]] = data.frame(Outlier_selected)


   writexl::write_xlsx(Outlier_DF, "Outlier_DF.xlsx")

   # Save list of outliers as .csv as well, so that it is accepted by DESEQ_NORM

   write.csv(data.frame(Outlier_selected), "Outliers_Selected.csv", row.names=FALSE)


   setwd(Output_file_path)

   return()

}













