# DEG per time point and Pathway Enrichment. All comparisons are calculated by this function. The output folder structure creates a folder for all molecules
DEG_FUNCTION = function(Output_file_path, List_contrasts_Path, DEG_Method = 'DESeq',MH_Method = 'BH', AlphaHC = 0.1,  padjval = 0.2 ,LogFoldThrs_VolPlot = 1, AdjustDeSeq =c("Plate.id")){
    require(ggplot2)
    require(ggrepel)
    require(DESeq2)
    #require(ddpca)
    require(foreach)
    require(BioMark)
    require(stringr)
    require("openxlsx")
    require(dplyr)
    require(doParallel)
    require(BiocParallel)
    require(parallel)

    # source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/Functions_Invivo.R")
    # Intro message
    print("Getting DEG...")

    # Import data
    setwd(Output_file_path)
    Data_file = "./DESEQ_NORM_QCNORM/DESeq_Norm.RData"
    load(Data_file)

    print("Only contrasts allowed after removing outliers will be calculated")
    List_contrasts = read.xlsx(xlsxFile = List_contrasts_Path)
    List_contrasts = List_contrasts %>% filter(treat %in% Metadata$CoarseCondition & untreat %in% Metadata$CoarseCondition)

    drugsInContrasts = unique(sapply(List_contrasts$treat, function(x){str_split(x,'_')[[1]][1]}))
    drugs = unique(Metadata$Treatment)
    drugsInContrasts = drugsInContrasts[drugsInContrasts %in% drugs]

    # Create Folder
    Name_folder =  "DEG_FUNCTION"
    dir.create(Name_folder)
    setwd(Name_folder)



    # DEG per molecule. The .RData data frame contains all genes with Log2Fold changes, the .xlsx have been filtered by padjval
    for(drug in drugsInContrasts){

        Name_folder_mol =drug
        dir.create(Name_folder_mol)
        setwd(Name_folder_mol)
        print(drug)


        List_contrasts_mol =  List_contrasts %>% filter(treat %in% c(List_contrasts$treat[grep(drug, List_contrasts$treat)]))

        DiffGenes = list()
        #wb <- createWorkbook()

        f.treat_ttest = function(ri,curData,curConc, ControlRef) {
            Results = list()
            gene = rownames(curData)[ri]
            y = curData[ri,]
            x = curConc
            df = data.frame(y,x)
            if(var(y) != 0){
                Ttest = t.test(y ~ x, data = df, var.equal = FALSE)
                p.value = Ttest$p.value
            }else{
                p.value = 1

            }

            log2FoldChange = mean(y[x != ControlRef]) - mean(y[x == ControlRef])

            Results[['ModelInfo']] = c(genes = gene, p.value = p.value, log2FoldChange = log2FoldChange)
            return(Results)
        }
        f.treat_negbin = function(ri,curData,curConc, ControlRef) {
            Results = list()
            gene = rownames(curData)[ri]
            y = curData[ri,]
            x = curConc
            df = data.frame(y,x)
            if(var(y) != 0){
                fit = glm(y ~ x, data = df, family = poisson(link = "log"))
                log2FoldChange = coef(summary(fit))['xVehicle','Estimate']
                p.value = coef(summary(fit))['xVehicle' , 'Pr(>|z|)']
            }else{
                p.value = 1
                log2FoldChange = 0

            }



            Results[['ModelInfo']] = c(genes = gene, p.value = p.value, log2FoldChange = log2FoldChange)
            return(Results)
        }
        f.treat_lm = function(ri,curData,curConc, ControlRef) {
            Results = list()
            gene = rownames(curData)[ri]
            y = log2(curData[ri,]+1)
            x = curConc
            df = data.frame(y,x)
            if(var(y) != 0){
                fit = lm(y ~ x, data = df)
                log2FoldChange = coef(summary(fit))['xVehicle','Estimate']
                p.value = coef(summary(fit))['xVehicle' , 'Pr(>|t|)']
            }else{
                p.value = 1
                log2FoldChange = 0

            }



            Results[['ModelInfo']] = c(genes = gene, p.value = p.value, log2FoldChange = log2FoldChange)
            return(Results)
        }
        f.treat_Wilcox = function(ri,curData,curConc, ControlRef) {
            Results = list()
            gene = rownames(curData)[ri]
            y = curData[ri,]
            x = curConc
            df = data.frame(y,x)
            if(var(y) != 0){
                test = wilcox.test(y ~ x, data = df)
                log2FoldChange = test$statistic
                p.value = test$p.value
            }else{
                p.value = 1
                log2FoldChange = 0

            }



            Results[['ModelInfo']] = c(genes = gene, p.value = p.value, log2FoldChange = log2FoldChange)
            return(Results)
        }

        # List for normalized matrices for DESeq_Cons method
        Norm_List = list()
        for(i in 1:nrow(List_contrasts_mol)){
            print(List_contrasts_mol$treat[i])
            Name_contr = paste(List_contrasts_mol$treat[i],'VS',List_contrasts_mol$untreat[i],sep = '_')
            #Name_contr = paste(List_contrasts_mol$treat[i])
            print(Name_contr)
            # Filter necessary data for comparison from Metadata
            Metadata_temp = Metadata %>% dplyr::filter(CoarseCondition %in% c(List_contrasts_mol$treat[i],List_contrasts_mol$untreat[i]))
            NormCounts_temp = NormCounts[,match(Metadata_temp$Sample_name, colnames(NormCounts))]

            curConc = Metadata_temp$CoarseCondition
            ControlRef = List_contrasts_mol$untreat[i]

            if(DEG_Method == 'T-Test'){

                print('T-test method for DEG')
                treatRes = foreach (ri = 1:nrow(NormCounts_temp),.export=c('f.treat_ttest'),
                                    .packages=c()) %dopar% {
                                        res <- f.treat_ttest(ri,curData = NormCounts_temp,curConc, ControlRef)
                                        return(res)

                                    }
                ModelInfo = lapply(treatRes, function(x){return(x$ModelInfo)})
                ModelInfo = as.data.frame(do.call(rbind,ModelInfo))
                rownames(ModelInfo) = ModelInfo$genes
                ModelInfo$p.value = as.numeric(ModelInfo$p.value)

            }else if(DEG_Method == 'DESeq'){

                print('DESeq method')
                res = results(deseqObj,contrast = c("CoarseCondition",List_contrasts_mol$treat[i], List_contrasts_mol$untreat[i] ) )
                ModelInfo = data.frame(res)
                ModelInfo$genes = rownames(ModelInfo)
            }else if(DEG_Method == 'DESeq_Cons'){

                print('DESeq_Cons method')
                # Download counts
                Count = counts(deseqObj)
                Count_temp = Count[,match(Metadata_temp$Sample_name, colnames(Count))]
                Count_temp = as.data.frame(Count_temp)
                Design = "~ CoarseCondition"
                # To be modified if more than one factor is desired to be adjusted for
                if(AdjustDeSeq != 1 ){

                        Design = paste(Design, AdjustDeSeq, sep = ' + ' )



                }

                dds_temp <- DESeqDataSetFromMatrix(countData = Count_temp, colData = Metadata_temp, design = as.formula(Design) )
                startTime <- Sys.time()
                deseqObj_temp = DESeq(dds_temp,
                                 parallel = TRUE,
                                 fitType = "parametric",
                                 BPPARAM=MulticoreParam(detectCores() - 4))

                print("Using Vst ")
                NormCounts_temp = getVarianceStabilizedData(deseqObj_temp)
                Norm_List[[Name_contr]] = NormCounts_temp
                res = results(deseqObj_temp,contrast = c("CoarseCondition",List_contrasts_mol$treat[i], List_contrasts_mol$untreat[i] ) )
                ModelInfo = data.frame(res)
                ModelInfo$genes = rownames(ModelInfo)

            }else{

                stop('DEG method not valid')
            }


            if( MH_Method == 'High_Cr'){

                print('High criticism method for multiple hypothesis testing')
                # Mix both because it can get misleading results only HCthresh when all the HCi are negative (which would mean veery null or something).
                # First check whether HC detection is rejected or not and then estimate the significant pvalues.

                HCdet = HCdetection(ModelInfo$p.value, alpha = AlphaHC, pvalcut = NA)
                if(HCdet$H == 0){
                    print('HC did not reject null hypothesis')
                    ModelInfo$pvaj = 0

                }else if(HCdet$H == 1){
                    print('HC rejecting null hypothesis')
                    Ind = HCthresh(ModelInfo$p.value, alpha = AlphaHC, plotit = TRUE)
                    ModelInfo$pvaj =0
                    ModelInfo$pvaj[Ind] = 1
                }
                ModelInfo_sig = ModelInfo %>% dplyr::filter(pvajSig == 1)

            }else if(MH_Method == 'BH'){
                if(DEG_Method == 'T-Test'){
                    ModelInfo$padj =p.adjust(as.numeric(ModelInfo$p.value), method='BH')
                    ModelInfo_sig = ModelInfo %>% dplyr::filter(padj < padjval)
                }else if(DEG_Method %in%  c('DESeq', 'DESeq_Cons', 'limma_voom')){

                    ModelInfo_sig = ModelInfo %>% filter(padj < padjval)
                    print(nrow(ModelInfo_sig))
                }

                # Plot volcano plot
                VP = Volcano_plot_padj(res_Vol = ModelInfo, padj_thr = padjval, logFC_thr = LogFoldThrs_VolPlot, Title = Name_contr , maxOver = 20,xLabel = 'Log2Fold',yLabel = '-Log10(P-adj)')

                png(filename=paste("VolcanoPlot", Name_contr,".png", sep = '_'), width = 1000, height = 1000, res=100)
                print(VP)
                dev.off()

            }else{

                stop("Select a valid Multiple hypothesis method")
            }



            DiffGenes[[Name_contr]] = ModelInfo
            # save .xlsx
            write.xlsx(ModelInfo, paste('Cln', unique(Metadata$Cell_line), 'Mol',drug, 'comp', Name_contr ,'.xlsx', sep = '_'))

            #addWorksheet(wb, Name_contr)
            #writeData(wb, sheet = Name_contr, x = ModelInfo_sig)


        }

        if(DEG_Method == 'DESeq_Cons'){

            save(Norm_List,file = paste('NormCounts', drug,'.RData' , sep = '_') )
        }

        save(DiffGenes, file = paste('DEG', drug,'.RData' , sep = '_'))
        #saveWorkbook(wb, paste('DEG', drug,'.xlsx', sep = '_'), overwrite = TRUE)

        setwd(paste(Output_file_path,Name_folder, sep = '/'))


    }
    setwd(Output_file_path)
    return()
}

# DEG per condition and the folder structure is per comparison. Volcano plots are being generated in the report itself. Only one comparison is calculated by this function.
DEG_FUNCTION_DA = function(comp_vect, saveDir,Metadata, count_data, DEG_Method = 'DESeq_Cons', MH_Method = 'BH', AdjustDeSeq =c("Plate.id"), padj_thr_gene = 0.2, logFC_thrs_gene = 0){

    require(DESeq2)
    require(dplyr)
    require(foreach)
    require(BioMark)
    require(stringr)
    require("openxlsx")
    require(doParallel)
    require(BiocParallel)
    require(parallel)
    require(edgeR)
    require(limma)

    # source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/Functions_Invivo.R")
    # Intro message
    print("Getting DEG...")


    # DEG per molecule. The .RData data frame contains all genes with Log2Fold changes, the .xlsx have been filtered by padjval

    f.treat_ttest = function(ri,curData,curConc, ControlRef) {
        Results = list()
        gene = rownames(curData)[ri]
        y = curData[ri,]
        x = curConc
        df = data.frame(y,x)
        if(var(y) != 0){
            Ttest = t.test(y ~ x, data = df, var.equal = FALSE)
            p.value = Ttest$p.value
        }else{
            p.value = 1

        }

        log2FoldChange = mean(y[x != ControlRef]) - mean(y[x == ControlRef])

        Results[['ModelInfo']] = c(genes = gene, p.value = p.value, log2FoldChange = log2FoldChange)
        return(Results)
    }
    f.treat_negbin = function(ri,curData,curConc, ControlRef) {
        Results = list()
        gene = rownames(curData)[ri]
        y = curData[ri,]
        x = curConc
        df = data.frame(y,x)
        if(var(y) != 0){
            fit = glm(y ~ x, data = df, family = poisson(link = "log"))
            log2FoldChange = coef(summary(fit))['xVehicle','Estimate']
            p.value = coef(summary(fit))['xVehicle' , 'Pr(>|z|)']
        }else{
            p.value = 1
            log2FoldChange = 0

        }



        Results[['ModelInfo']] = c(genes = gene, p.value = p.value, log2FoldChange = log2FoldChange)
        return(Results)
    }
    f.treat_lm = function(ri,curData,curConc, ControlRef) {
        Results = list()
        gene = rownames(curData)[ri]
        y = log2(curData[ri,]+1)
        x = curConc
        df = data.frame(y,x)
        if(var(y) != 0){
            fit = lm(y ~ x, data = df)
            log2FoldChange = coef(summary(fit))['xVehicle','Estimate']
            p.value = coef(summary(fit))['xVehicle' , 'Pr(>|t|)']
        }else{
            p.value = 1
            log2FoldChange = 0

        }



        Results[['ModelInfo']] = c(genes = gene, p.value = p.value, log2FoldChange = log2FoldChange)
        return(Results)
    }
    f.treat_Wilcox = function(ri,curData,curConc, ControlRef) {
        Results = list()
        gene = rownames(curData)[ri]
        y = curData[ri,]
        x = curConc
        df = data.frame(y,x)
        if(var(y) != 0){
            test = wilcox.test(y ~ x, data = df)
            log2FoldChange = test$statistic
            p.value = test$p.value
        }else{
            p.value = 1
            log2FoldChange = 0

        }



        Results[['ModelInfo']] = c(genes = gene, p.value = p.value, log2FoldChange = log2FoldChange)
        return(Results)
    }

    #print(comp_vect)

    Contrast_name = paste0(comp_vect['treat'],"_vs_",comp_vect['untreat'])
    # Create dir
    dir.create(file.path(saveDir, "contrasts", Contrast_name), recursive = TRUE)
    resfilename <- file.path(saveDir, "contrasts", Contrast_name, "DEG_df.RDS")

    # Filter necessary data for comparison from Metadata
    Metadata$CoarseCondition = gsub("-", "_",Metadata$CoarseCondition)
    Metadata_temp = Metadata %>% dplyr::filter(CoarseCondition %in% c(comp_vect['treat'],comp_vect['untreat']))
    NormCounts_temp = NormCounts[,match(Metadata_temp$Sample_name, colnames(NormCounts))]
    count_data_temp =count_data[,match(Metadata_temp$Sample_name, colnames(count_data))]

    curConc = Metadata_temp$CoarseCondition
    ControlRef = comp_vect['untreat']

    # T-test needs to be adapted to new pipeline
    if(DEG_Method == 'T-Test'){

            print('T-test method for DEG')
            treatRes = foreach (ri = 1:nrow(NormCounts_temp),.export=c('f.treat_ttest'),
                                .packages=c()) %dopar% {
                                    res <- f.treat_ttest(ri,curData = NormCounts_temp,curConc, ControlRef)
                                    return(res)

                                }
            ModelInfo = lapply(treatRes, function(x){return(x$ModelInfo)})
            ModelInfo = as.data.frame(do.call(rbind,ModelInfo))
            rownames(ModelInfo) = ModelInfo$genes
            ModelInfo$p.value = as.numeric(ModelInfo$p.value)

    }else if(DEG_Method == 'DESeq_Cons'){

            #print('DESeq_Cons method')
            deseq_filename <- file.path(saveDir, "contrasts", Contrast_name, "DEseqObj_cons.RDS")
            norm_filename <- file.path(saveDir, "contrasts", Contrast_name, "Norm_cons.RDS")

            # Download counts

            Count_temp = as.data.frame(count_data_temp)
            Design = "~ CoarseCondition"
            # To be modified if more than one factor is desired to be adjusted for
            if(AdjustDeSeq != "" ){
                if(length(unique(Metadata_temp[[AdjustDeSeq]])) > 1){
                Design = paste(Design, AdjustDeSeq, sep = ' + ' )
                }
            }


            dds_temp <- DESeqDataSetFromMatrix(countData = Count_temp, colData = Metadata_temp, design = as.formula(Design) )
            deseqObj_cons = DESeq(dds_temp,
                                  parallel = TRUE,
                                  fitType = "parametric",
                                  BPPARAM=MulticoreParam(detectCores() - 4))

            #print("Using Vst ")
            NormCounts_cons = getVarianceStabilizedData(deseqObj_cons)
            # Preparing for heatmap

            res = results(deseqObj_cons,contrast = c("CoarseCondition",comp_vect['treat'], comp_vect['untreat'] ), parallel = TRUE)
            # shrink the lfcs
            #res <-  lfcShrink(deseqObj_cons, res = res, type = "ashr")
            res <- tibble(symbol = rownames(res),
                          log2FC = res$log2FoldChange,
                          pvalue = res$pvalue,
                          padj = res$padj,
                          lfcSE = res$lfcSE,
                          baseMean=res$baseMean)
            res<- res %>% arrange(padj)


            saveRDS(res, resfilename)
            saveRDS(deseqObj_cons, deseq_filename)
            saveRDS(NormCounts_cons, norm_filename)


        }else if(DEG_Method == 'limma_voom'){

            print('limma-voom method')
            voom_filename <- file.path(saveDir, "contrasts", Contrast_name, "voomObj_cons.RDS")
            norm_filename <- file.path(saveDir, "contrasts", Contrast_name, "Norm_cons.RDS")

            if(AdjustDeSeq != "" ){
                stop("limma_voom method has not been adapted to adjust for batch")
            }

            # To be modified if more than one factor is desired to be adjusted for

            design <- model.matrix(~ 0 + CoarseCondition, data=Metadata_temp)
            colnames(design) = gsub("CoarseCondition", "",colnames(design)  )
            Count_temp = as.data.frame(count_data_temp)
            keep <- rowSums(Count_temp) > 1
            Count_temp <- Count_temp[keep,]
            kept_genes <- rownames(Count_temp)[keep]
            dge <- DGEList(counts = Count_temp )
            dge <- calcNormFactors(dge, method="TMM")

            v <- voom(dge, design, plot = TRUE)
            fit <- lmFit(v, design)
            # Make the contrast correctly
            contr_expr <- paste0(comp_vect['treat'], " - ", comp_vect['untreat'])
            contr <- makeContrasts(contrasts = contr_expr, levels = colnames(coef(fit)))
            tmp <- contrasts.fit(fit, contr)
            tmp <- eBayes(tmp)
            res0 <- topTable(tmp, number=Inf)

            res <- tibble(symbol = rownames(res0),
                         log2FC = res0$logFC,
                         pvalue = res0$P.Value,
                         padj = res0$adj.P.Val,
                         AveExpr = res0$AveExpr)

            # Get normalized counts
            NormCounts_cons = v$E

            saveRDS(res, resfilename)
            saveRDS(fit, voom_filename)
            saveRDS(NormCounts_cons, norm_filename)



        }else{

            stop('DEG method not valid')
        }



        # if( MH_Method == 'High_Cr'){
        #
        #     print('High criticism method for multiple hypothesis testing')
        #     # Mix both because it can get misleading results only HCthresh when all the HCi are negative (which would mean veery null or something).
        #     # First check whether HC detection is rejected or not and then estimate the significant pvalues.
        #
        #     HCdet = HCdetection(ModelInfo$p.value, alpha = AlphaHC, pvalcut = NA)
        #     if(HCdet$H == 0){
        #         print('HC did not reject null hypothesis')
        #         ModelInfo$pvaj = 0
        #
        #     }else if(HCdet$H == 1){
        #         print('HC rejecting null hypothesis')
        #         Ind = HCthresh(ModelInfo$p.value, alpha = AlphaHC, plotit = TRUE)
        #         ModelInfo$pvaj =0
        #         ModelInfo$pvaj[Ind] = 1
        #     }
        #     ModelInfo_sig = ModelInfo %>% filter(pvajSig == 1)
        #
        # }else if(MH_Method == 'BH'){
        #     if(DEG_Method == 'T-Test'){
        #         ModelInfo$padj =p.adjust(as.numeric(ModelInfo$p.value), method='BH')
        #         ModelInfo_sig = ModelInfo %>% filter(padj < padjval)
        #     }else if(DEG_Method %in%  c('DESeq', 'DESeq_Cons')){
        #
        #         ModelInfo_sig = ModelInfo %>% filter(padj < padjval)
        #         print(nrow(ModelInfo_sig))
        #     }
        #
        #
        #
        # }else{
        #
        #     stop("Select a valid Multiple hypothesis method")
        # }

        # # Plot heatmaps
        # Filenames
        Top_up_filename <- file.path(saveDir, "contrasts", Contrast_name, "Top_up_DEGs.RDS")
        Top_up_filename_jpeg <- file.path(saveDir, "contrasts", Contrast_name, "Top_up_DEGs.jpeg")
        Top_dn_filename <- file.path(saveDir, "contrasts", Contrast_name, "Top_dn_DEGs.RDS")
        Top_dn_filename_jpeg <- file.path(saveDir, "contrasts", Contrast_name, "Top_dn_DEGs.jpeg")
        All_DEG_filename_pdf <- file.path(saveDir, "contrasts", Contrast_name, "All_DEGs.pdf")

        mat_anno <- Metadata_temp %>% arrange(Treatment) %>% dplyr::select(Treatment, Stimulant_used, Treatment_conc, Outliers, Sample_name) %>% column_to_rownames(var = "Sample_name")
        Heat_dat = NormCounts_cons[, match(rownames(mat_anno), colnames(NormCounts_cons))]

        # Top upregulated

        #anno_colors <- list(sample = setNames(gg_color_hue(length(unique(Metadata_temp$sample))), unique(Metadata_temp$sample)))
        top_genes <- res %>% dplyr::filter(log2FC > logFC_thrs_gene & padj < padj_thr_gene) %>% top_n(-50, padj) %>% pull(symbol)

        if(length(top_genes) > 0){
            Heat_dat_filtered <- Heat_dat[match(top_genes, rownames(Heat_dat)), , drop = FALSE]

            lenchar = round(max(nchar(colnames(Heat_dat_filtered)))/10)

            if(length(top_genes)>1){
                require(dendsort)
                sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), isReverse = TRUE))
                cluster_rows <- sort_hclust(hclust(dist( Heat_dat_filtered)))
                Hup = pheatmap(Heat_dat_filtered,
                         color = colorRampPalette(c("blue", "white", "red"))(200),
                         show_rownames = TRUE,
                         annotation_col = mat_anno,
                         border_color = NA,
                         fontsize = 10,
                         scale = "row",
                         cluster_cols = T,
                         cluster_rows = cluster_rows,
                         fontsize_row = 8)
             }else {
                Hup = pheatmap(Heat_dat_filtered,
                         color = colorRampPalette(c("blue", "white", "red"))(200),
                         show_rownames = TRUE,
                         annotation = mat_anno,
                         border_color = NA,
                         fontsize = 10,
                         scale = "row",
                         cluster_cols = T,
                         cluster_rows = F,
                         fontsize_row = 8)
             }

            jpeg(file = Top_up_filename_jpeg, width = 10, height = (lenchar + nrow(Heat_dat_filtered) * 0.2), units = "in", res = 300)
            print(Hup)
            dev.off()

            saveRDS(Hup, Top_up_filename)


        }



        # # Top downregulated
        top_genes_dn <- res %>% dplyr::filter(log2FC < -logFC_thrs_gene & padj < padj_thr_gene) %>% top_n(-50, padj) %>% pull(symbol)

        if(length(top_genes_dn) > 0){

            Heat_dat_filtered <- Heat_dat[match(top_genes_dn, rownames(Heat_dat)), , drop = FALSE]
            lenchar = round(max(nchar(colnames(Heat_dat_filtered)))/10)
            if(length(top_genes_dn)>1){
                require(dendsort)
                sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), isReverse = TRUE))
                cluster_rows <- sort_hclust(hclust(dist( Heat_dat_filtered)))
                Hdn = pheatmap(Heat_dat_filtered,
                               color = colorRampPalette(c("blue", "white", "red"))(200),
                               annotation = mat_anno,
                               show_rownames = TRUE,
                               border_color = NA,
                               fontsize = 10,
                               scale = "row",
                               cluster_cols = T,
                               cluster_rows = cluster_rows,
                               fontsize_row = 8)
            } else {
                Hdn = pheatmap(Heat_dat_filtered,
                               color = colorRampPalette(c("blue", "white", "red"))(200),
                               show_rownames = TRUE,
                               annotation = mat_anno,
                               border_color = NA,
                               fontsize = 10,
                               scale = "row",
                               cluster_cols = T,
                               cluster_rows = F,
                               fontsize_row = 8)
            }

            jpeg(file = Top_dn_filename_jpeg, width = 10, height = (lenchar + nrow(Heat_dat_filtered) * 0.2), units = "in", res = 300)
            print(Hdn)
            dev.off()

            saveRDS(Hdn, Top_dn_filename)

        }



        # All DEGs
        deg_genes <- res %>% dplyr::filter(padj < padj_thr_gene) %>% pull(symbol)
        if(length(deg_genes) > 0){
            Heat_dat_filtered <- Heat_dat[match(deg_genes, rownames(Heat_dat)), , drop = FALSE]
            if(length(deg_genes)>1){
                require(dendsort)
                sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), isReverse = TRUE))
                cluster_rows <- sort_hclust(hclust(dist( Heat_dat_filtered)))
                Hall = pheatmap(Heat_dat_filtered,
                               color = colorRampPalette(c("blue", "white", "red"))(200),
                               annotation = mat_anno,
                               show_rownames = TRUE,
                               border_color = NA,
                               fontsize = 10,
                               scale = "row",
                               cluster_cols = T,
                               cluster_rows = cluster_rows,
                               fontsize_row = 8)
            } else {
                Hall = pheatmap(Heat_dat_filtered,
                               color = colorRampPalette(c("blue", "white", "red"))(200),
                               show_rownames = T,
                               annotation = mat_anno,
                               border_color = NA,
                               fontsize = 10,
                               scale = "row",
                               cluster_cols = T,
                               cluster_rows = F,
                               fontsize_row = 8)
            }

            #jpeg(file = All_DEG_filename_jpeg, width = 12, height = (ceiling(length(deg_genes)/10)+5), units = "in", res = 300)
            pdf(file = All_DEG_filename_pdf, width = 12, height = (ceiling(length(deg_genes)/10)+5))
            print( Hall)
            dev.off()

        }





  return(res)
}
