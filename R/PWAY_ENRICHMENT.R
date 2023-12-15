# DEG per time point and Pathway Enrichment
PWAY_ENRICHMENT = function(Output_file_path, WITH_REVERSE = TRUE, REVERSE_DRUG = 'Sham', REVERSE_TIME = "72", padjval = 0.2, LogFoldThrs = 1,Pway_qvalThrs = 0.2, ORGANISM = 'Mouse', Target_list_path, DirData){
    library(ggvenn)
    library(ggplot2)
    library(ggrepel)
    library("AnnotationDbi")
    library("org.Mm.eg.db")
    library("clusterProfiler")
    library("pathview")
    library(DESeq2)
    library(eulerr)
    library("org.Hs.eg.db")

    # Function to add Chemoproteomic targets present within an enriched pathway
    Add_Chemo_Enrich = function(drug, Target_list, enrichedPathways){
        # This function adds a column in the enrichedPathways data frame with the chemoproteomics targets associated to that molecule present in the pathway

        # drug: Treatment for which enrichment was generated
        # Target_list: Target list from chemoproteomics
        # enrichedPathways: data frame with enriched pathways from clusterProfiler

        # Add column for candidate targets
        Target_mol = Target_list %>% dplyr::filter(drugs == drug)
        Target = Target_mol$entrezid
        names(Target) = Target_mol$gene

        WithTar = c()
        # Data frame V1
        if(!is.null(Target)){
            for(pway in enrichedPathways$ID){

                GenesPway = GeneSet %>% dplyr::filter(from == pway)
                WithTar_t = paste(names(Target[Target %in% GenesPway$to]), collapse = ', ')
                WithTar = c(WithTar, WithTar_t)
            }
            enrichedPathways$Target = WithTar
        }else{

            enrichedPathways$Target = "na"
        }

        return(enrichedPathways)



    }


    # Checking parameters if FILTER_REVERSE = FALSE does not matter what REVERSE_CONTROL is, cause it is not used. If FILTER_REVERSE = TRUE, the length of REVERSE_CONTROL
    # should be bigger than zero.
    # Data bases for pathway enrichment
    if(ORGANISM == 'Mouse'){
        # Mouse
        KEGG_GeneSet_path = paste(DirData,"KEGG_All_mouse_230922.RData", sep = "/")
        REACTOME_GeneSet_path = paste(DirData,"REACTOME_All_mouse_230924.RData", sep = "/")
        Species = 'mmu'
    }else if(ORGANISM == 'Human'){
        KEGG_GeneSet_path = paste(DirData,"KEGG_All_human_231208.RData", sep = "/")
        REACTOME_GeneSet_path = paste(DirData,"REACTOME_All_human_231209.RData", sep = "/")
        Species = 'hsa'

    }else{

        stop("select a valid organism")
    }

    # Download list of targets from Chemoproteomics
    if(file.exists(Target_list_path)){
        load(Target_list_path)
    }else{
        stop("Provide a valid path for chemoproteomics list")
    }

    # Create Folder
    Name_folder =  "PWAY_ENRICHMENT"
    dir.create(Name_folder)
    setwd(Name_folder)


    # Find drugs with DEG
    drugs = list.files("../DEG_FUNCTION")



    if(WITH_REVERSE == TRUE){

        if(!REVERSE_DRUG %in% drugs){

            stop("REVERSE_DRUG does not match any drug in the experiment. REVERSE_DRUG is needed for WITH_REVERSE = TRUE")
        }else{

            FILTERS = c('TRUE', 'FALSE')

        }


    }else{


        FILTERS = c('FALSE')
    }


    for(FILTER_REVERSE in FILTERS){


            if(FILTER_REVERSE == TRUE){
                # Create Folder
                Name_folder_filter =  "PWAY_REVERSE_TRUE"
                dir.create(Name_folder_filter)
                setwd(Name_folder_filter)

                if(FILTER_REVERSE == TRUE){

                    drugs = drugs[drugs != REVERSE_DRUG]
                }
            }else if (FILTER_REVERSE == FALSE){

                # Create Folder
                Name_folder_filter =  "PWAY_REVERSE_FALSE"
                dir.create(Name_folder_filter)
                setwd(Name_folder_filter)


            }


            for(drug in drugs){

                Name_folder_mol = drug
                dir.create(Name_folder_mol)
                setwd(Name_folder_mol)
                print(drug)

                if(FILTER_REVERSE == TRUE){
                    # Download Sham DEG
                    NameFile = paste("DEG", REVERSE_DRUG, '.RData', sep = '_')
                    DEG_file_path = paste("../../../DEG_FUNCTION",REVERSE_DRUG,NameFile, sep = '/')
                    load(DEG_file_path)
                    DiffGenes_healthy = DiffGenes
                }

                # Download drug DEG
                NameFile = paste("DEG", drug, '.RData', sep = '_')
                DEG_file_path = paste("../../../DEG_FUNCTION",drug,NameFile, sep = '/')
                load(DEG_file_path)


                # List of comparisons or pathway enrichment
                Listdrug = names(DiffGenes)[grep(drug,names(DiffGenes))]

                DiffGenesPways_kegg = list()
                DiffGenesPways_Reactome = list()
                DiffGenesPways_GO = list()
                EnrichedPathwaysDf_kegg = data.frame("ID" = character(), "Description" = character(), "GeneRatio" = character, "BgRatio" = character,
                                                     "pvalue" = numeric(), "p.adjust" = numeric(), "qvalue" = numeric(), "geneID" = character(),
                                                     "Count" = integer(), "Target" = character(), "comp" = character() )

                EnrichedPathwaysDf_Reactome = data.frame("ID" = character(), "Description" = character(), "GeneRatio" = character, "BgRatio" = character,
                                                         "pvalue" = numeric(), "p.adjust" = numeric(), "qvalue" = numeric(), "geneID" = character(),
                                                         "Count" = integer(), "Target" = character(), "comp" = character())


                EnrichedPathwaysDf_GO = data.frame("ID" = character(), "Description" = character(), "GeneRatio" = character, "BgRatio" = character,
                                                   "pvalue" = numeric(), "p.adjust" = numeric(), "qvalue" = numeric(), "geneID" = character(),
                                                   "Count" = integer(), "Target" = character(), "comp" = character() )






                for(comp in Listdrug){

                    print(comp)

                    if(FILTER_REVERSE == 'TRUE'){

                            # Plot Venn diagram of genes reversing disease

                            Df_mol = DiffGenes[[comp]]
                            Df_mol_filt = Df_mol %>% filter(padj < padjval & abs(log2FoldChange) > LogFoldThrs)

                            if(nrow(Df_mol_filt) == 0){

                                print(paste(comp, "does not have DEG"))
                                next
                            }

                            if(!is.null(REVERSE_TIME)){

                                Healthy_Name = paste(REVERSE_DRUG, REVERSE_TIME, sep = '_')
                                Healthy_Comp = paste(Healthy_Name,'_VS',str_split(comp, 'VS')[[1]][2],sep = '')
                                if(!Healthy_Comp %in% names(DiffGenes_healthy)){

                                    print(paste(Healthy_Comp, "was not calculated in DEG step", sep = ' '))
                                    next
                                }else{

                                    Df_mol_sham = DiffGenes_healthy[[Healthy_Comp]]
                                    Df_mol_sham_filt = Df_mol_sham %>% filter(padj < padjval & abs(log2FoldChange) > LogFoldThrs)

                                }


                            }else{

                                stop("Under construction")
                            }


                            Name1 = str_split(comp, 'VS')[[1]][1]
                            Name2 = str_split(Healthy_Comp, 'VS')[[1]][1]
                            x = list(Df_mol_filt$genes, Df_mol_sham_filt$genes)
                            names(x) = c(Name1, Name2)
                            if(drug != REVERSE_DRUG){
                                P = plot(euler(x), quantities = TRUE)
                                png(filename=paste("VennDiagram", comp,'Filterlogfold', LogFoldThrs,".png", sep = '_'), width = 1000, height = 800, res=100)
                                print(P)
                                dev.off()
                            }
                            # Build data frame

                            Df_mol_filt_sel = Df_mol_filt %>% dplyr::select(genes, log2FoldChange)
                            colnames(Df_mol_filt_sel) = c('genes', "log2F_mol")
                            Df_mol_sham_filt_sel = Df_mol_sham_filt %>% dplyr::select(genes, log2FoldChange)
                            colnames(Df_mol_sham_filt_sel) = c('genes', "log2F_sham")

                            Inter = merge(Df_mol_filt_sel,Df_mol_sham_filt_sel, by = 'genes', all = FALSE)
                            Inter$ProdSign = sign(Inter$log2F_mol) * sign(Inter$log2F_sham)

                            InterSumm = Inter %>% group_by(ProdSign) %>% summarise(count=n())
                            InterSumm$Sim = 'Same'
                            InterSumm$Sim[InterSumm$ProdSign == -1] = 'Opposite'

                            p<-ggplot(data=InterSumm, aes(x=Sim, y=count)) +
                                geom_bar(stat="identity")+
                                theme_minimal()+
                                ggtitle(paste("Intersection ", Name1 ,' and ',Name2 , sep = ''))

                            png(filename=paste("OverlapGenes_dir", comp,'Filterlogfold', LogFoldThrs,".png", sep = '_'), width = 1000, height = 800, res=100)
                            print(p)
                            dev.off()

                            Name_folder_L2 = comp
                            dir.create(Name_folder_L2)
                            setwd(Name_folder_L2)

                            Df_mol_filt$Log2Sham = Df_mol_sham_filt$log2FoldChange[match(rownames(Df_mol_filt), rownames(Df_mol_sham_filt))]
                            Df_mol_filt = Df_mol_filt[!is.na(Df_mol_filt$Log2Sham),]
                            Df_mol_filt$Simil = sign(Df_mol_filt$log2FoldChange) * sign(Df_mol_filt$Log2Sham)
                            Df_mol_filt = Df_mol_filt %>% filter(Simil ==1)



                    }else{

                        Df_mol = DiffGenes[[comp]]
                        Df_mol_filt = Df_mol %>% filter(padj < padjval & abs(log2FoldChange) > LogFoldThrs)

                        if(nrow(Df_mol_filt) == 0){

                            print(paste(comp, "does not have DEG"))
                            next
                        }

                        Name_folder_L2 = comp
                        dir.create(Name_folder_L2)
                        setwd(Name_folder_L2)


                    }


                    # Add entrezid to the data frame
                    if(ORGANISM == 'Mouse'){

                        Df_mol_filt$entrez <- mapIds(org.Mm.eg.db, keys=rownames(Df_mol_filt), column="ENTREZID", keytype="SYMBOL", multiVals="first")


                    }else if(ORGANISM == 'Human'){

                        Df_mol_filt$entrez <- mapIds(org.Hs.eg.db, keys=rownames(Df_mol_filt), column="ENTREZID", keytype="SYMBOL", multiVals="first")

                    }

                    naIdx = which(is.na(Df_mol_filt$entrez))
                    if (length(naIdx) > 0) {
                        Df_mol_filt = Df_mol_filt[-naIdx,]
                    }

                    # If number of gene
                    if(nrow(Df_mol_filt)< 3){
                        print('Not enough DEG genes for pathway enrichment after converting to symbols')
                        setwd('..')
                        next
                    }


                    #KEGG

                    Name_folder_L3 = "KEGG"
                    dir.create(Name_folder_L3)
                    setwd(Name_folder_L3)


                    print('Enrich KEGG')

                    load(KEGG_GeneSet_path)
                    try(res <- enricher(Df_mol_filt$entrez, TERM2GENE = GeneSet,TERM2NAME=GeneSet_Desc))


                    if(!is.null(res)){
                        if(nrow(summary(res)) >1 & max(summary(res)$Count) > 1){
                            enrichedPathways = res@result %>% filter(p.adjust < Pway_qvalThrs)

                            # Fixing notation
                            resultTemp = res@result
                            if(ORGANISM == 'Mouse'){

                               resultTemp$Description = sapply(resultTemp$Description, function(x){str_split(x,' - Mus')[[1]][1]})
                               res@result = resultTemp
                               enrichedPathways$Description = sapply(enrichedPathways$Description, function(x){str_split(x,' - Mus')[[1]][1]})

                               # This data frame for the cnetplot
                               resx <- setReadable(res, "org.Mm.eg.db", 'ENTREZID')

                            }

                            # Add chemoproteomics candidates to the pathways enriched in the data frame
                            enrichedPathways = Add_Chemo_Enrich(drug, Target_list, enrichedPathways)
                            # Save data frame of enriched pathways
                            enrichedPathways$comp = comp
                            EnrichedPathwaysDf_kegg = rbind(EnrichedPathwaysDf_kegg,enrichedPathways )
                            save(enrichedPathways, file = paste("EnrichPathways","Pway_qvalThrs", Pway_qvalThrs,".RData", sep = "_"))


                            # Dotplot
                            png(file="DotPlot_Enrich.png",width = 800, height = 800, res = 100)
                            print(try(dotplot(res, showCategory=20, font.size = 12) + ggtitle(paste(drug, comp, sep = ' '))))
                            dev.off()

                            # cnetplot
                            # It works on entrezID
                            genelist = Df_mol_filt$log2FoldChange
                            names(genelist) = Df_mol_filt$entrez
                            png(file="CnetPlot_Enrich.png",width = 800, height = 800, res = 100)
                            if(ORGANISM == "Mouse"){
                               print(try(cnetplot(resx, foldChange=genelist, showCategory = 10)))
                            }else{

                                print(try(cnetplot(res, foldChange=genelist, showCategory = 10)))
                            }
                            dev.off()

                            # Plot Kegg diagrams

                            #Problematic_Pways = c("mmu04723", "mmu00512", "mmu01212")
                            Problematic_Pways = c()
                            print("Got till here!")
                            for (ti in 1:length(enrichedPathways$ID)){

                                curKegg = enrichedPathways$ID[ti]
                                print(curKegg)
                                curKegg_Desc = enrichedPathways$Description[ti]
                                gids = unlist(str_split(enrichedPathways$geneID[ti], '/'))
                                GeneFold = Df_mol_filt$log2FoldChange[which(Df_mol_filt$entrez %in% gids)]
                                names(GeneFold) = Df_mol_filt$entrez[which(Df_mol_filt$entrez %in% gids)]
                                DiffGenesPways_kegg[[paste(comp,curKegg_Desc,sep='_' )]] = Df_mol_filt %>% filter(entrez %in% gids) %>% dplyr::select("log2FoldChange", "pvalue", "padj", "genes", "entrez")
                                print(GeneFold)
                                if (!(enrichedPathways$ID[ti] %in% Problematic_Pways)){
                                    try(pathview(gene.data  = GeneFold,
                                                 pathway.id = curKegg,
                                                 limit      = list(gene=max(abs(GeneFold)), cpd=1),
                                                 kegg.dir = '.',
                                                 kegg.native= TRUE,
                                                 out.suffix = paste(drug,curKegg_Desc,sep = "_"),
                                                 species    = Species,
                                                 res = 100,
                                                 low=list(gene="steelblue"),
                                                 high=list(gene="aquamarine3")))
                                }
                            }



                        }else{
                            print("KEGG enrichment is not very reliable, therefore not saved")
                        }
                    }
                    setwd("..")


                    #REACTOME

                    Name_folder_L3 = "REACTOME"
                    dir.create(Name_folder_L3)
                    setwd(Name_folder_L3)


                    print('Enrich REACTOME')
                    # Get Reactome from mouse

                    load(REACTOME_GeneSet_path)
                    try(res <- enricher(Df_mol_filt$entrez, TERM2GENE = GeneSet,TERM2NAME=GeneSet_Desc))

                    if(!is.null(res)){
                      if(nrow(summary(res)) >1 & max(summary(res)$Count) > 1){
                        enrichedPathways = res@result %>% filter(qvalue < Pway_qvalThrs)

                        # Add chemoproteomics candidates to the pathways enriched in the data frame
                        enrichedPathways = Add_Chemo_Enrich(drug, Target_list, enrichedPathways)

                        # Save enriched pathways
                        enrichedPathways$comp = comp
                        EnrichedPathwaysDf_Reactome = rbind(EnrichedPathwaysDf_Reactome,enrichedPathways )
                        save(enrichedPathways, file = paste("EnrichPathways","Pway_qvalThrs", Pway_qvalThrs,".RData", sep = "_"))



                        # Plots
                        Dot = dotplot(res, showCategory=20, font.size = 12) + ggtitle(paste(drug, comp, sep = ' '))
                        # Dotplot
                        png(file="DotPlot_Enrich.png",width = 800, height = 800, res = 100)
                        print(Dot)
                        dev.off()

                        # cnetplot
                        # Adding symbol
                        if(ORGANISM == 'Mouse'){
                            # This data frame for the cnetplot
                            resx <- setReadable(res, "org.Mm.eg.db", 'ENTREZID')

                        }
                        genelist = Df_mol_filt$log2FoldChange
                        names(genelist) = Df_mol_filt$entrez

                        png(file="CnetPlot_Enrich.png",width = 800, height = 800, res = 100)
                        if(ORGANISM == "Mouse"){
                            print(try(cnetplot(resx, foldChange=genelist, showCategory = 10)))
                        }else{

                            print(try(cnetplot(res, foldChange=genelist, showCategory = 10)))
                        }
                        dev.off()

                        Problematic_Pways = c()
                        for (ti in 1:length(enrichedPathways$ID)){
                            if (!(enrichedPathways$ID[ti] %in% Problematic_Pways)){
                                curKegg = enrichedPathways$ID[ti]
                                print(curKegg)
                                curKegg_Desc = enrichedPathways$Description[ti]
                                gids = unlist(str_split(enrichedPathways$geneID[ti], '/'))
                                GeneFold = Df_mol_filt$log2FoldChange[which(Df_mol_filt$entrez %in% gids)]
                                names(GeneFold) = Df_mol_filt$entrez[which(Df_mol_filt$entrez %in% gids)]

                                DiffGenesPways_Reactome[[paste(comp,curKegg_Desc,sep='_' )]] = Df_mol_filt %>% filter(entrez %in% gids) %>% dplyr::select("log2FoldChange", "pvalue", "padj", "genes", "entrez")

                            }
                        }

                      }else{

                          print("REACTOME enrichment is not very reliable, therefore not saved")

                       }
                    }

                    setwd("..")



                    # GO Enrichment

                    Name_folder_L3 = "GO"
                    dir.create(Name_folder_L3)
                    setwd(Name_folder_L3)

                    print('Enrich GO')

                    if(ORGANISM == 'Mouse'){

                        try(res <- enrichGO(Df_mol_filt$entrez, OrgDb = "org.Mm.eg.db", ont="ALL",readable=TRUE, qvalueCutoff = Pway_qvalThrs))

                    }else if(ORGANISM == 'Human'){


                        try(res <- enrichGO(Df_mol_filt$entrez, OrgDb = 'org.Hs.eg.db', ont="ALL",readable=TRUE, qvalueCutoff = Pway_qvalThrs))
                    }else{

                        stop("Select a valid organism")
                    }

                    if(!is.null(res)){
                       if(nrow(summary(res)) >1 & max(summary(res)$Count) > 1){

                        enrichedPathways = res@result %>% filter(p.adjust < Pway_qvalThrs)

                        # Add chemoproteomics candidates to the pathways enriched in the data frame
                        enrichedPathways = Add_Chemo_Enrich(drug, Target_list, enrichedPathways)

                        # Save enriched pathways
                        enrichedPathways$comp = comp
                        EnrichedPathwaysDf_GO = rbind(EnrichedPathwaysDf_GO,enrichedPathways )
                        save(enrichedPathways, file = paste("EnrichPathways","Pway_qvalThrs", Pway_qvalThrs,".RData", sep = "_"))


                        # Plots

                        # Dotplot
                        png(file="DotPlot.png", width = 800, height = 800, res = 100)
                        try(print(dotplot(res, showCategory=20, font.size = 10, split = "ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")))
                        dev.off()

                        # cnetplot
                        # Adding symbol
                        if(ORGANISM == 'Mouse'){
                            # This data frame for the cnetplot
                            resx <- setReadable(res, "org.Mm.eg.db", 'ENTREZID')

                        }
                        genelist = Df_mol_filt$log2FoldChange
                        names(genelist) = Df_mol_filt$entrez

                        png(file="CnetPlot_Enrich.png",width = 800, height = 800, res = 100)
                        if(ORGANISM == "Mouse"){
                            print(try(cnetplot(resx, foldChange=genelist, showCategory = 10)))
                        }else{

                            print(try(cnetplot(res, foldChange=genelist, showCategory = 10)))
                        }
                        dev.off()

                        Problematic_Pways = c()
                        for (ti in 1:length(enrichedPathways$ID)){
                            if (!(enrichedPathways$ID[ti] %in% Problematic_Pways)){
                                curKegg = enrichedPathways$ID[ti]
                                print(curKegg)
                                curKegg_Desc = enrichedPathways$Description[ti]
                                gids = unlist(str_split(enrichedPathways$geneID[ti], '/'))
                                GeneFold = Df_mol_filt$log2FoldChange[which(Df_mol_filt$Gene %in% gids)]
                                names(GeneFold) = Df_mol_filt$entrez[which(Df_mol_filt$Gene %in% gids)]
                                print(GeneFold)
                                DiffGenesPways_GO[[paste(comp,curKegg_Desc,sep='_' )]] = Df_mol_filt %>% filter(genes %in% gids) %>% dplyr::select("log2FoldChange", "pvalue", "padj", "genes", "entrez")

                            }
                        }
                       }else{

                           print("GO enrichment is not very reliable, therefore not saved")

                      }
                    }
                    setwd("..")
                    setwd("..")

                }


                save(EnrichedPathwaysDf_kegg,EnrichedPathwaysDf_Reactome,EnrichedPathwaysDf_GO, DiffGenesPways_kegg, DiffGenesPways_Reactome, DiffGenesPways_GO ,
                     file = paste('EnrichPathways_AllGeneSets', 'FILTER_REVERSE', FILTER_REVERSE,'.RData', sep = '_'))
                setwd('..')

            }

            setwd('..')
    }

}








