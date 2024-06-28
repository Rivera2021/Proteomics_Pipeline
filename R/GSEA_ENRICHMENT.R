GSEA_ENRICHMENT = function(DirDataForPipeline, Output_file_path, QCNORM, ORGANISM, pval, List_contrasts_Path   ){

    library(knitr)
    library(ggplot2)
    library(ggpubr)
    library(rstatix)
    library(tidyverse)
    library(ggpubr)
    library(openxlsx)
    library(org.Hs.eg.db)
    library(DESeq2)
    library(BiocParallel)
    library(parallel)
    library(pheatmap)
    library(fgsea)
    library(stringr)
    library(dplyr)

    # Function for finding gsea enrich pathways
    # NOTE: The size within the gsea enriched data frame (I think) is the number of genes in the geneset that are actually contained in the universe of our data set.
    # NOTE2:
    Enrich_gsea = function(gene_list, Gmt_file, pvalue) {
        set.seed(54321)
        library(dplyr)
        library(fgsea)


        if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
            warning("We had to sort genes")
            gene_list = sort(gene_list, decreasing = TRUE)
        }
        myGS = fgsea::gmtPathways(Gmt_file)

        gsea_Res <- fgsea::fgsea(pathways = myGS,
                                 stats = gene_list,
                                 minSize=10, ## minimum size of gene set
                                 maxSize=400, ## maximum size of gene set
                                 nperm=10000) %>%
            as.data.frame() %>%
            dplyr::filter(padj < !! pvalue) %>%
            arrange(desc(NES))
        message(paste("Number of enriched pathways =", nrow(gsea_Res)))

        # Collapse pathways
        concise_pathways = collapsePathways(data.table::as.data.table(gsea_Res),
                                            pathways = myGS,
                                            stats = gene_list)
        gsea_Res = gsea_Res[gsea_Res$pathway %in% concise_pathways$mainPathways, ]
        message(paste("Number of enriched pathways after collapse =", nrow(gsea_Res)))

        gsea_Res$Enrichment = ifelse(gsea_Res$NES > 0, "Up regulated", "Down regulated")
        gsea_Res = gsea_Res[order(gsea_Res$padj, decreasing = FALSE),]

        filtRes = rbind(head(gsea_Res, n = 20))

        total_up = sum(gsea_Res$Enrichment == "Up regulated")
        total_down = sum(gsea_Res$Enrichment == "Down regulated")
        header = paste0("Top 10: Up=", total_up,", Down=",    total_down, ")")

        colos = setNames(c("firebrick2", "dodgerblue2"),
                         c("Up regulated", "Down regulated"))

        g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
            geom_point( aes(fill = Enrichment, size = size), shape=21) +
            scale_fill_manual(values = colos ) +
            scale_size_continuous(range = c(2,10)) +
            geom_hline(yintercept = 0) +
            coord_flip() +
            labs(x="Pathway", y="Normalized Enrichment Score",
                 title=header)

        output = list("Results" = gsea_Res, "Plot" = g1)
        return(output)
    }

    plot_geneset_clusters = function( Gsea_Res, Gmt_file, min.sz = 3, Title){
        library(ggplot2)
        library(ggrepel)
        library(stringr)
        library(dynamicTreeCut)

        myGS = fgsea::gmtPathways(Gmt_file)
        df = matrix(nrow=nrow(Gsea_Res), ncol = nrow(Gsea_Res), data = 0)
        rownames(df) = colnames(df) = Gsea_Res$pathway

        for ( i in 1:nrow(Gsea_Res)) {
            genesI =  unlist(myGS[names(myGS) == Gsea_Res$pathway[i] ])
            for (j in 1:nrow(Gsea_Res)) {
                genesJ = unlist(myGS[names(myGS) == Gsea_Res$pathway[j] ])
                ## Jaccards distance  1 - (intersection / union )
                overlap = sum(!is.na(match(genesI, genesJ )))
                jaccards = overlap / length(unique(c(genesI, genesJ) ))
                df[i,j] = 1-jaccards
            }
        }

        ## Cluster nodes using dynamic tree cut, for colors
        distMat = as.dist(df)
        dendro = hclust(distMat, method = "average" )
        clust = dynamicTreeCut::cutreeDynamicTree( dendro, minModuleSize = min.sz )
        ## Note: in dynamicTreeCut, cluster 0, is a garbage cluster for things that dont cluster, so we remove it

        Gsea_Res$Cluster = clust
        #Gsea_Res = Gsea_Res[Gsea_Res$Cluster != 0, ]

        ## select gene sets to label for each clusters
        relevant = Gsea_Res %>% group_by( Cluster ) %>% top_n(wt = abs(size), n = 1) %>% .$pathway
        # determine cluster order for plotting
        clust_ords = Gsea_Res %>% group_by( Cluster ) %>% summarise("Average" = NES ) %>% arrange(desc(Average)) %>% .$Cluster %>% unique

        Gsea_Res$Cluster = factor(Gsea_Res$Cluster, levels = clust_ords)

        Gsea_Res$Label = ""
        Gsea_Res$Label[Gsea_Res$pathway %in% relevant ] = Gsea_Res$pathway[Gsea_Res$pathway %in% relevant ]
        #Gsea_Res$Label = str_remove(Gsea_Res$Label, "HALLMARK_")
        Gsea_Res$Label = tolower(Gsea_Res$Label)

        g1 = ggplot(Gsea_Res, aes(x = Cluster, y = NES, label = Label )) +
            geom_jitter( aes(color = Cluster,  size = size), alpha = 0.8, height = 0, width = 0.2 ) +
            scale_size_continuous(range = c(0.5,5)) +
            geom_text_repel( force = 2, max.overlaps = Inf) +
            ggtitle(Title)

        Results = list(Plt = g1, relevant = relevant)

        return(Results)
    }



    Gene_set_names = c('HALLMARK')
    if(ORGANISM == 'Mouse'){

        stop('GSEA for mouse under construction for Mouse')

    }else if(ORGANISM == 'Human'){

        Hallmark_gmt_file = "h.all.v2023.2.Hs.symbols.gmt"
        Hallmark_GeneSet_gmt_path = paste(DirDataForPipeline, Hallmark_gmt_file , sep = "/")
        Species = 'hsa'

    }else{

        stop("select a valid organism")
    }


    # Download normalized matrix

    if(QCNORM == "PRE_QCNORM"){

        Data_Norm = "./DESEQ_NORM/DESeq_Norm.RData"
        load(Data_Norm)


    }else if(QCNORM == "POST_QCNORM"){

        Data_Norm = "./DESEQ_NORM_QCNORM/DESeq_Norm.RData"
        load(Data_Norm)
        # Take out samples that won't pair because an entire control condition was removed
        #Metadata = Metadata %>% filter(!timeColl %in% TimeToRemove)
        #NormCounts = NormCounts[,match(Metadata$Sample_name, colnames(NormCounts))]


    }


    # Download list of constrasts
    List_contrasts = read.xlsx(xlsxFile = List_contrasts_Path)

    # Create Folder
    setwd(Output_file_path)

    Name_folder =  "GSEA_ENRICHMENT"
    dir.create(Name_folder)
    setwd(Name_folder)


    # Find treatments DEG

    if(QCNORM == "PRE_QCNORM"){

        # Find drugs with DEG
        drugs = list.files("../DEG_FUNCTION")



    }else if(QCNORM == "POST_QCNORM"){

        # Find drugs with DEG
        drugs = list.files("../DEG_FUNCTION_QCNORM")

    }

    # # GSEA per Time per dose
    #
    # ControlRef = 'DMSO'
    #
    # MOLEC = unique(Metadata$Treatment)
    # #MOLEC = MOLEC[!MOLEC %in% c('DMSO', 'NONE')]

    for(mol in drugs){

        Name_folder_mol = mol
        dir.create(Name_folder_mol)
        setwd(Name_folder_mol)
        print(mol)


        if(QCNORM == "PRE_QCNORM"){

            NameFile = paste("DEG", mol, ".RData", sep = '_')
            DEG_file_path = paste("../../DEG_FUNCTION",mol,NameFile, sep = '/')
            load(DEG_file_path)



        }else if(QCNORM == "POST_QCNORM"){

            NameFile = paste("DEG", mol, ".RData", sep = '_')
            DEG_file_path = paste("../../DEG_FUNCTION_QCNORM",mol,NameFile, sep = '/')
            load(DEG_file_path)

        }

        MapDf = data.frame(NamesDiff = names(DiffGenes))
        MapDf$Mol = sapply(MapDf$NamesDiff, function(x){str_split(x, '_')[[1]][1]})
        MapDf$Dose = sapply(MapDf$NamesDiff, function(x){str_split(x, '_')[[1]][2]})
        MapDf$Time = sapply(MapDf$NamesDiff, function(x){str_split(x, '_')[[1]][3]})
        MapDf$Stim = sapply(MapDf$NamesDiff, function(x){str_split(x, '_')[[1]][4]})

        for(time in unique(MapDf$Time)){


            MapDf_time = MapDf %>% filter(Time == time)
            DiffGenes_time = DiffGenes[MapDf_time$NamesDiff]

            for(stim in unique(MapDf_time$Stim)){

                MapDf_time_stim = MapDf_time %>% filter(Stim == stim)
                DiffGenes_time_stim = DiffGenes_time[MapDf_time_stim$NamesDiff]

                MapDf_time_stim = MapDf_time_stim[order(as.numeric(MapDf_time_stim$Dose)),]

                for(dose in unique(MapDf_time_stim$Dose)){


                    MapDf_time_stim_dose = MapDf_time_stim %>% filter(Dose == dose)
                    DiffGenes_time_stim_dose = DiffGenes_time_stim[MapDf_time_stim_dose$NamesDiff][[1]]
                    DiffGenes_time_stim_dose$gsea_in = sign(DiffGenes_time_stim_dose$log2FoldChange) * (-log10(DiffGenes_time_stim_dose$pvalue))
                    DiffGenes_time_stim_dose = DiffGenes_time_stim_dose[order(DiffGenes_time_stim_dose$gsea_in, decreasing = TRUE),]

                    Name_folder_mol_comp = MapDf_time_stim_dose$NamesDiff
                    dir.create(Name_folder_mol_comp)
                    setwd(Name_folder_mol_comp)
                    print(Name_folder_mol_comp)


                    gene_list  = DiffGenes_time_stim_dose$gsea_in
                    names(gene_list) = rownames(DiffGenes_time_stim_dose)
                    gene_list = gene_list[!is.na(gene_list)]

                    # HALLMARK
                    Gmt_file = Hallmark_GeneSet_gmt_path
                    Hall_gsea = Enrich_gsea(gene_list, Gmt_file, pval)
                    GSEA_df = Hall_gsea$Results
                    Plt = Hall_gsea$Plot
                    #save(GSEA_df, file = paste(MapDf_time_stim_dose$NamesDiff, ".RData", sep = '_'))
                    ggsave(filename=paste("DotPLot",MapDf_time_stim_dose$NamesDiff, ".pdf", sep = '_'), plot= Plt, width = 8,height = 6, units = 'in')
                    ggsave(filename=paste("DotPLot",MapDf_time_stim_dose$NamesDiff, ".jpeg", sep = '_'), plot=Plt, width = 8,height = 6, units = 'in')

                    # PLOT CLUSTERS OF PATHWAYS

                    Res_up = GSEA_df[GSEA_df$NES > 0, ]
                    Res_down = GSEA_df[GSEA_df$NES < 0, ]
                    Relev = c()
                    if(nrow(Res_up) > 1){
                        Resclus_up = plot_geneset_clusters( Gsea_Res = Res_up, Gmt_file, min.sz = 3, Title = 'GSEA clusters: Up-regulated')
                        Relev_up = Resclus_up$relevant
                        ggsave(filename=paste("Cluster_Up",MapDf_time_stim_dose$NamesDiff, ".jpeg", sep = '_'), plot = Resclus_up$Plt, width = 8,height = 6, units = 'in')
                        Relev = c(Relev,Relev_up)
                    }

                    if(nrow(Res_down) > 1){
                        Resclus_down = plot_geneset_clusters( Gsea_Res = Res_down , Gmt_file, min.sz = 3, Title = 'GSEA clusters: Down-regulated')
                        Relev_down = Resclus_down$relevant
                        ggsave(filename=paste("Cluster_Down",MapDf_time_stim_dose$NamesDiff, ".jpeg", sep = '_'), plot = Resclus_down$Plt, width = 8,height = 6, units = 'in')
                        Relev = c(Relev,Relev_down)
                    }
                    save(GSEA_df,Relev, file = paste(MapDf_time_stim_dose$NamesDiff, ".RData", sep = '_'))

                    # PLOT GENE HEATMAPS OF GENES INVOLVED IN GSEA PATHWAYS
                    Comp_Prep = List_contrasts[List_contrasts$treat == MapDf_time_stim_dose$NamesDiff, ]
                    for(pway in GSEA_df$pathway){

                        GeneRelev = GSEA_df$leadingEdge[GSEA_df$pathway == pway]
                        Metadata_comp = Metadata %>% filter(CoarseCondition %in% c(Comp_Prep$treat, Comp_Prep$untreat) )
                        if(mol %in% c('DMSO', 'NONE')){

                            Metadata_comp = Metadata_comp[order(Metadata_comp$Stimulant_used), ]
                            annotation_col = Metadata_comp %>% dplyr::select(Stimulant_used, Outliers)
                            rownames(annotation_col) = Metadata_comp$Sample_name


                        }else{

                            Metadata_comp = Metadata_comp[order(Metadata_comp$Treatment), ]
                            annotation_col = Metadata_comp %>% dplyr::select(Treatment, Outliers)
                            rownames(annotation_col) = Metadata_comp$Sample_name
                        }



                        NormCounts_comp = NormCounts[match(GeneRelev[[1]], rownames(NormCounts)), match(Metadata_comp$Sample_name, colnames(NormCounts)),drop=FALSE]
                        NormCounts_comp_zscore = t(sapply(as.data.frame(t(NormCounts_comp)), function(x) {(x-mean(x))/sd(x)}))
                        colnames(NormCounts_comp_zscore) = colnames(NormCounts_comp)


                        jpeg(file = paste("Heatmap", MapDf_time_stim_dose$NamesDiff,pway,'.jpeg',sep ='_' ), width = 6, height = 6, units = "in", res = 300)

                         P = pheatmap(as.matrix(NormCounts_comp_zscore), cluster_cols=FALSE, cluster_rows=FALSE, annotation_col = annotation_col, fontsize_row = 4, fontsize_col = 4, angle_col = 270)
                         print(P)

                        dev.off()


                    }

                    setwd('..')



                }



            }


        }

        setwd('..')




    }


    setwd(Output_file_path)




}
