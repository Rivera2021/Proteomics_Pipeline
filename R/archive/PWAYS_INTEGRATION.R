PWAYS_INTEGRATION = function(Output_file_path,Target_list_path, ORGANISM = 'Mouse', DirPipeline_Data, Control ){

    library(openxlsx)
    library(dplyr)
    library(stringr)

    # Functions
    Is_Tar_in_Pway = function(TargetCand,Target_list, EnrichedPways, GeneSet){
        # This function adds a column in the enrichedPathways data frame with the chemoproteomics targets associated to that molecule present in the pathway
        # TargetCand: The target we want to evaluate presence in pathway
        # Target_list: Target list from chemoproteomics
        # enrichedPathways: data frame with enriched pathways from clusterProfiler
        # GeneSet: KEGG Gene set

        # Add column for candidate targets
        Target = Target_list$entrezid[match(TargetCand, Target_list$gene)]
        names(Target) = TargetCand

        WithTar = c()

        for(pway in EnrichedPways$ID){

            GenesPway = GeneSet %>% dplyr::filter(from == pway)
            WithTar_t = paste(names(Target[Target %in% GenesPway$to]),collapse = '')
            WithTar = c(WithTar, WithTar_t)
        }
        EnrichedPways$WithTarget = WithTar


        return(EnrichedPways)



    }

    # Create Folder
    setwd(Output_file_path)
    Name_folder =  "PWAYS_INTEGRATION"
    dir.create(Name_folder)
    setwd(Name_folder)


    # CREATE SUMMARY SPREAD SHEET FOR DEG

    Mols = list.files(paste('../DEG_FUNCTION', sep = '/'))
    # Main Folder
    MainFolder_DEG = paste('../DEG_FUNCTION', sep = '/')

    Diff_Genes_mol = list()
    Diff_Genes_Tot = list()
    for(mol in Mols){
        print(mol)
        Filename = paste("DEG", mol,".RData", sep = '_')
        load(paste(MainFolder_DEG, mol, Filename, sep = '/'))
        print(names(DiffGenes))
        DiffGenes = DiffGenes[grep(mol, names(DiffGenes))]
        DiffGenes_filt = list()
        for(i in 1: length(DiffGenes)){
            print(i)
            df = DiffGenes[[i]]
            #df$Gene = rownames(df)
            #df = df %>% filter(padj <0.2) %>% dplyr::select(Gene,log2FoldChange)
            df = df %>% filter(padj <0.2) %>% dplyr::select(genes,log2FoldChange)
            colnames(df) = c("Gene",names(DiffGenes)[i] )
            DiffGenes_filt[[names(DiffGenes)[i]]] = df

        }
        Diff_Genes_mol[[mol]] = DiffGenes_filt %>% purrr::reduce(full_join, by = "Gene")

    }


    Diff_Genes_Tot = Diff_Genes_mol %>% purrr::reduce(full_join, by = "Gene")

    wb <- createWorkbook()
    addWorksheet(wb, "DEG_All_Comparisons")
    writeData(wb, sheet = 'DEG_All_Comparisons', x = Diff_Genes_Tot)

    saveWorkbook(wb, "DEG_All_Comparisons.xlsx", overwrite = TRUE)


   rm(Diff_Genes_Tot,Diff_Genes_mol,DiffGenes_filt, DiffGenes)





    # CREATE SUMMARY SPREAD SHEETS FOR ENRICHMENT
    FolderType = list.files("../PWAY_ENRICHMENT")
    Gene_set_names = c('KEGG', 'REACTOME', 'GO')
    Ntop_GO = 30

    for(Fold in FolderType){

        Mols = list.files(paste('../PWAY_ENRICHMENT',Fold, sep = '/'))

            for(Gene_set_name in Gene_set_names){

                EnrichedPathwaysDf_All = data.frame("ID" = character(),"Description" = character(), "GeneRatio" = character(),   "BgRatio"=character() , "pvalue" = numeric() ,"p.adjust" = numeric(),  "qvalue" = numeric(), "geneID" = character() , "Count" = integer() ,"Target"= character() ,  "comp" = character() ,"mol" = character())


                for(mol in Mols){

                    Filestemp = list.files(paste('../PWAY_ENRICHMENT',Fold,mol, sep = '/'))
                    EnrichPwayFile = Filestemp[grep("EnrichPathways_AllGeneSets",Filestemp)]
                    load(paste('../PWAY_ENRICHMENT',Fold,mol,EnrichPwayFile, sep = '/'))

                    if(Gene_set_name == 'KEGG'){

                        EnrichedPathwaysDf_KEGG$mol = mol
                        EnrichedPathwaysDf_All = rbind(EnrichedPathwaysDf_All, EnrichedPathwaysDf_KEGG)


                    }else if(Gene_set_name == 'REACTOME'){

                        EnrichedPathwaysDf_REACTOME$mol = mol
                        EnrichedPathwaysDf_All = rbind(EnrichedPathwaysDf_All, EnrichedPathwaysDf_REACTOME)


                    }else if(Gene_set_name == 'GO' ){

                        EnrichedPathwaysDf_GO$mol = mol
                        EnrichedPathwaysDf_All = rbind(EnrichedPathwaysDf_All, EnrichedPathwaysDf_GO)

                    }



                }

                # Save .RData frame
                save(EnrichedPathwaysDf_All, file = paste("EnrichedPways_All", Gene_set_name,Fold,'.RData',sep = '_') )

                # Add extra filtering for GO
                if(Gene_set_name == 'GO'){
                    # Rebuild Go data frame with only the top 25 pathways
                    EnrichedPathwaysDf_All = data.frame("ID" = character(),"Description" = character(), "GeneRatio" = character(),   "BgRatio"=character() , "pvalue" = numeric() ,"p.adjust" = numeric(),  "qvalue" = numeric(), "geneID" = character() , "Count" = integer() ,"Target"= character() ,  "comp" = character() ,"mol" = character())
                    for(mol in Mols){

                        Filestemp = list.files(paste('../PWAY_ENRICHMENT',Fold,mol, sep = '/'))
                        EnrichPwayFile = Filestemp[grep("EnrichPathways_AllGeneSets",Filestemp)]
                        load(paste('../PWAY_ENRICHMENT',Fold,mol,EnrichPwayFile, sep = '/'))


                        EnrichedPathwaysDf_GO_filt  = EnrichedPathwaysDf_GO %>% group_by(comp) %>% arrange(desc(Count)) %>% slice_head(n = Ntop_GO)
                        EnrichedPathwaysDf_GO_filt$mol = mol
                        EnrichedPathwaysDf_All = rbind(EnrichedPathwaysDf_All, EnrichedPathwaysDf_GO_filt)





                    }

                }




                # Remove columns not needed in excel files

                EnrichedPathwaysDf_All = subset(EnrichedPathwaysDf_All, select = -c(qvalue,NameLongRef, NameRef) )

                # # EnrichedPathwaysDf_All$Namesheet = paste('Genes',EnrichedPathwaysDf_All$CompTemp, gsub('/','_',EnrichedPathwaysDf_All$geneID), '.png',sep = '_')
                # EnrichedPathwaysDf_All$Namesheet = sapply(EnrichedPathwaysDf_All$Namesheet, function(x){substr(x,1,28)})
                # seq <- 1:nrow(EnrichedPathwaysDf_All)
                # seq = seq %% 17

                # Shortened NameRefUniq to fit the restrictions of excel name

                EnrichedPathwaysDf_All$NameRefUniq2 = sapply(EnrichedPathwaysDf_All$NameRefUniq, function(x){gsub(paste(Control,'_', sep = ''),'',x)})
                EnrichedPathwaysDf_All$NameRefUniq2 = gsub(':','_',EnrichedPathwaysDf_All$NameRefUniq2)


                # Create excel workbook
                wb <- createWorkbook()
                Mainsheet = paste("EnrichedPways", Gene_set_name,sep = '_')
                addWorksheet(wb, Mainsheet)
                writeData(wb, sheet = Mainsheet, x = EnrichedPathwaysDf_All)


                # Add sheets for gene expression. Unique combinations of genes

                for(g in unique(EnrichedPathwaysDf_All$NameRefUniq)){

                    MOL = EnrichedPathwaysDf_All$mol[match(g, EnrichedPathwaysDf_All$NameRefUniq)]
                    COMP = EnrichedPathwaysDf_All$comp[match(g, EnrichedPathwaysDf_All$NameRefUniq)]
                    Namesheet = EnrichedPathwaysDf_All$NameRefUniq2[match(g, EnrichedPathwaysDf_All$NameRefUniq)]

                    # Add gene expression plot

                    addWorksheet(wb, Namesheet)

                    Fileplot = paste('Genes',g,'.png', sep = '_')
                    Fileplot_path = paste(Output_file_path, 'PWAY_ENRICHMENT', Fold,MOL, COMP,Gene_set_name, Fileplot, sep = '/')

                    if(Gene_set_name %in% c('REACTOME', 'KEGG', 'GO')){
                        if(file.exists(Fileplot_path)){
                            insertImage(wb,sheet =Namesheet, file = Fileplot_path,
                                        width = 8,
                                        height = 6,
                                        startRow = 3,
                                        startCol = 1,
                                        units = "in",
                                        dpi = 500
                            )
                        }
                    }
                    # Add link to main page
                    writeFormula(wb, Namesheet,
                                 startRow = 1, startCol = 15,
                                 x = makeHyperlinkString(sheet = Mainsheet, row = 1, col = 1, text = "Back to main")
                    )

                    # Add data frame with descriptions

                    Filestemp = list.files(paste('../PWAY_ENRICHMENT',Fold,MOL, sep = '/'))
                    EnrichPwayFile = Filestemp[grep("EnrichPathways_AllGeneSets",Filestemp)]
                    load(paste('../PWAY_ENRICHMENT',Fold,MOL,EnrichPwayFile, sep = '/'))

                    FileDesc = paste('Genes',g, sep = '_')
                    #DiffUniqueGenesPwaysDf_Mol_KEGG[[FileDesc]]

                    if(Gene_set_name == 'KEGG'){
                        #print(FileDesc %in% names(DiffUniqueGenesPwaysDf_Mol_KEGG))
                        writeData(wb, sheet =  Namesheet, x = DiffUniqueGenesPwaysDf_Mol_KEGG[[FileDesc]],rowNames = FALSE,startRow = 3,
                                  startCol = 14)



                    }else if(Gene_set_name == 'REACTOME'){

                        #print(FileDesc %in% names(DiffUniqueGenesPwaysDf_Mol_REACTOME))
                        writeData(wb, sheet =  Namesheet, x = DiffUniqueGenesPwaysDf_Mol_REACTOME[[FileDesc]],rowNames = FALSE,startRow = 3,
                                  startCol = 14)


                    }else if(Gene_set_name == 'GO'){

                        print(FileDesc %in% names(DiffUniqueGenesPwaysDf_Mol_GO))
                        writeData(wb, sheet =  Namesheet, x = DiffUniqueGenesPwaysDf_Mol_GO[[FileDesc]],rowNames = FALSE,startRow = 3,
                                  startCol = 14)

                    }

                    # Add KEGG diagrams

                    Pway_desc = paste(MOL,EnrichedPathwaysDf_All$Description[match(g, EnrichedPathwaysDf_All$NameRefUniq)],sep = '_')
                    Folder_path = paste(Output_file_path, 'PWAY_ENRICHMENT', Fold,MOL, COMP,Gene_set_name, sep = '/')
                    FileDiagram = list.files(Folder_path)[grep(Pway_desc,list.files(Folder_path))]

                    FileDiagram_path = paste(Folder_path, FileDiagram, sep = '/')

                    if(Gene_set_name == 'KEGG'){
                        if(file.exists(FileDiagram_path)){
                            insertImage(wb,sheet =Namesheet, file = FileDiagram_path,
                                        width = 8,
                                        height = 6,
                                        startRow = 40,
                                        startCol = 1,
                                        units = "in",
                                        dpi = 500
                            )
                        }
                    }

                }

                # Add link from Main to plots
                row = 2
                for(nameid in EnrichedPathwaysDf_All$NameRefUniq2){


                    writeFormula(wb, Mainsheet,
                                 startRow = row, startCol = 13,
                                 x = makeHyperlinkString(sheet = nameid, row = 1, col = 1, text = nameid)
                    )

                    row = row + 1
                }


                saveWorkbook(wb, paste("EnrichedPways_All", Gene_set_name,Fold,'.xlsx',sep = '_'), overwrite = TRUE)
            }

    }




    # SUMMARY OF CHEMOPROTEOMICS
    print("Finding over-representation of chemoproteomics candidate targets")

    if(file.exists(Target_list_path)){
        load(Target_list_path)
    }else{
        stop("Provide a valid path for chemoproteomics list")
    }

    # Data bases for pathway enrichment
    if(ORGANISM == 'Mouse'){
        # Mouse
        KEGG_GeneSet_path = paste(DirPipeline_Data,"KEGG_All_mouse_230922.RData", sep = "/")
        REACTOME_GeneSet_path = paste(DirPipeline_Data,"REACTOME_All_mouse_230924.RData", sep = "/")
        GO_GeneSet_path = paste(DirPipeline_Data,"GO_All_mouse_240212.RData", sep = "/")

        Species = 'mmu'
    }else if(ORGANISM == 'Human'){
        KEGG_GeneSet_path = paste(DirPipeline_Data,"KEGG_All_human_231208.RData", sep = "/")
        REACTOME_GeneSet_path = paste(DirPipeline_Data,"REACTOME_All_human_231209.RData", sep = "/")
        GO_GeneSet_path = paste(DirPipeline_Data,"GO_All_human_240212.RData", sep = "/")
        Species = 'hsa'

    }else{

        stop("select a valid organism")
    }

    # Identify gene sets that contain a particular target
    TargetsTot = unique(Target_list$entrezid)

    # KEGG list
    load(KEGG_GeneSet_path)
    if(ORGANISM == 'Mouse'){

        GeneSet_Desc$to = sapply(GeneSet_Desc$to, function(x){str_split(x,' - Mus')[[1]][1]})

    }
    Chemo_Tar_KEGG = list()
    for(Tar in TargetsTot){

        Target_list_sel = Target_list %>% filter(entrezid == Tar)
        GeneSetWithTar = GeneSet %>% filter(to == Tar)
        GeneSetWithTar$Desc = GeneSet_Desc$to[match(GeneSetWithTar$from, GeneSet_Desc$from)]
        Chemo_Tar_KEGG[[unique(Target_list_sel$gene)[1]]] = GeneSetWithTar

    }

    # REACTOME list

    load(REACTOME_GeneSet_path)
    Chemo_Tar_REACTOME = list()
    for(Tar in TargetsTot){

        Target_list_sel = Target_list %>% filter(entrezid == Tar)
        GeneSetWithTar = GeneSet %>% filter(to == Tar)
        GeneSetWithTar$Desc = GeneSet_Desc$to[match(GeneSetWithTar$from, GeneSet_Desc$from)]
        Chemo_Tar_REACTOME[[unique(Target_list_sel$gene)[1]]] = GeneSetWithTar

    }


    # GO list

    load(GO_GeneSet_path)
    Chemo_Tar_GO = list()
    for(Tar in TargetsTot){

        Target_list_sel = Target_list %>% filter(entrezid == Tar)
        GeneSetWithTar = GeneSet %>% filter(to == Tar)
        GeneSetWithTar$Desc = GeneSet_Desc$to[match(GeneSetWithTar$from, GeneSet_Desc$from)]
        Chemo_Tar_GO[[unique(Target_list_sel$gene)[1]]] = GeneSetWithTar

    }


    # TABLE TO ESTIMATE OVER-REPRESENTATION OF A TARGET WITHIN ENRICHED PATHWAYS

    Gene_set_names = c('KEGG', 'REACTOME', 'GO')
    FolderType = list.files("../PWAY_ENRICHMENT")
    for(Fold in FolderType){

        drugs = list.files(paste("../PWAY_ENRICHMENT",Fold, sep = '/' ))



        for(Gene_set_name in Gene_set_names){

            print(Gene_set_name)

            # Create Data frame
            colnames_dyn = c()
            for(treat in drugs){
                colnames_t = c(paste(treat,'Tar_in', sep = '_'), paste(treat,'P_val', sep = '_'))
                colnames_dyn = c(colnames_dyn, colnames_t)

            }

            ColNames = c(c('CompCategory','Chemo_Target', 'Tar_in_GeneSet'), colnames_dyn)
            Chemo_Enrich = setNames(data.frame(matrix(ncol = length(ColNames), nrow = 0)), ColNames)
            Chemo_Enrich[] <- lapply(Chemo_Enrich, as.character)

            # Load enriched pathways big summary table
            load(paste("EnrichedPways_All", Gene_set_name, Fold,".RData", sep = '_') )

            # Add modification to table to make the task easier later
            EnrichedPathwaysDf_All$CompCategory = sapply(EnrichedPathwaysDf_All$comp, function(x){paste(str_split(x, '_')[[1]][2:length(str_split(x, '_')[[1]])], collapse = '_')})

            # Load KEGG Gene set

            if(Gene_set_name == 'KEGG'){
               load(KEGG_GeneSet_path)

            }else if(Gene_set_name == 'REACTOME'){

                load(REACTOME_GeneSet_path)
            }else if(Gene_set_name == 'GO'){

                load(GO_GeneSet_path)
            }

            Enriched_WithChemo =  EnrichedPathwaysDf_All %>% filter(Target != '')

            for(compo in unique(Enriched_WithChemo$CompCategory)){

                    Enriched_WithChemo_comp =  Enriched_WithChemo %>% filter(CompCategory == compo)
                    Targets_t = sapply(Enriched_WithChemo_comp$Target, function(x){str_split(x, ', ')[[1]] })
                    Targets_t = unique(unlist(Targets_t))
                    for(ChemoTar in Targets_t){


                        if(Gene_set_name == 'KEGG'){

                            Tar_in_KEGG = paste(as.character(nrow(Chemo_Tar_KEGG[[ChemoTar]])), '/', as.character(nrow(GeneSet_Desc)), sep = '')

                        }else if(Gene_set_name == 'REACTOME'){

                            Tar_in_REACTOME = paste(as.character(nrow(Chemo_Tar_REACTOME[[ChemoTar]])), '/', as.character(nrow(GeneSet_Desc)), sep = '')
                        }else if(Gene_set_name == 'GO'){

                            Tar_in_GO = paste(as.character(nrow(Chemo_Tar_GO[[ChemoTar]])), '/', as.character(nrow(GeneSet_Desc)), sep = '')
                        }


                        Tar_in_MOLS = c()
                        Pval_MOLS = c()
                        Joint_vect = c()
                        for(molec in drugs){

                            Enriched_comp_mol =  EnrichedPathwaysDf_All %>% filter(CompCategory == compo & mol == molec)

                            if(nrow(Enriched_comp_mol) > 0){

                                Enriched_comp_mol_WithTar = Is_Tar_in_Pway(TargetCand = ChemoTar,Target_list, EnrichedPways = Enriched_comp_mol, GeneSet)
                                Tar_in_Mol = nrow(Enriched_comp_mol_WithTar %>% filter(WithTarget != ''))

                            }else{

                                Tar_in_Mol = 0
                            }

                            No_Pways_Mol = nrow(Enriched_comp_mol)
                            Tar_in_MOLS[paste(molec,'Tar_in', sep = '_')] = paste(as.character(Tar_in_Mol), "/", as.character(No_Pways_Mol),sep = '')

                            # Hypergeometric distribution: Already sort of a cumulative if use phyper
                            x = Tar_in_Mol # P(X<= x)

                            if(Gene_set_name == 'KEGG'){

                               m = nrow(Chemo_Tar_KEGG[[ChemoTar]]) # Number of pathways containing the target

                            }else if(Gene_set_name == 'REACTOME'){

                                m = nrow(Chemo_Tar_REACTOME[[ChemoTar]])
                            }else if(Gene_set_name == 'GO'){

                                m = nrow(Chemo_Tar_GO[[ChemoTar]])
                            }

                            n = nrow(GeneSet_Desc) - m # Number of pathways not containing the target
                            k = No_Pways_Mol # Number of Pathways enriched
                            # x-1 cause I want the right tail from X >= x
                            Pval_MOLS[paste(molec,'P_val', sep = '_')] = sprintf("%e",  1 - phyper(x-1,m,n,k,lower.tail = TRUE ))
                            Joint_vect = c(Joint_vect,c(Tar_in_MOLS[paste(molec,'Tar_in', sep = '_')] ,  Pval_MOLS[paste(molec,'P_val', sep = '_')] ))



                        }

                        if(Gene_set_name == 'KEGG'){

                            Chemo_Enrich = Chemo_Enrich %>% add_row(data.frame(t(c('CompCategory' = compo,'Chemo_Target' = ChemoTar, 'Tar_in_GeneSet' = Tar_in_KEGG,Joint_vect))))
                            rm(m,n,k,x,Tar_in_KEGG)
                        }else if(Gene_set_name == 'REACTOME'){

                            Chemo_Enrich = Chemo_Enrich %>% add_row(data.frame(t(c('CompCategory' = compo,'Chemo_Target' = ChemoTar, 'Tar_in_GeneSet' = Tar_in_REACTOME,Joint_vect))))
                            rm(m,n,k,x,Tar_in_REACTOME)
                        }else if(Gene_set_name == 'GO'){

                            Chemo_Enrich = Chemo_Enrich %>% add_row(data.frame(t(c('CompCategory' = compo,'Chemo_Target' = ChemoTar, 'Tar_in_GeneSet' = Tar_in_GO,Joint_vect))))
                            rm(m,n,k,x,Tar_in_GO)
                        }


                    }
                }

            # Convert Pvalue columns in numeric:
            PvalCol = colnames(Chemo_Enrich)[grep('P_val', colnames(Chemo_Enrich))]
            for(colname in PvalCol){
                Chemo_Enrich[[colname]] = as.numeric(Chemo_Enrich[[colname]])
            }

            if(Gene_set_name == 'KEGG'){

                Chemo_Enrich_KEGG = Chemo_Enrich

            }else if(Gene_set_name == 'REACTOME'){


                Chemo_Enrich_REACTOME = Chemo_Enrich

            }else if(Gene_set_name == 'GO'){


                Chemo_Enrich_GO = Chemo_Enrich
            }


        }


        # Create excel workbook
        wb <- createWorkbook()
        addWorksheet(wb, "Chemo_Enrich_KEGG")
        writeData(wb, sheet = "Chemo_Enrich_KEGG", x = Chemo_Enrich_KEGG)
        addWorksheet(wb, "Chemo_Enrich_REACTOME")
        writeData(wb, sheet = "Chemo_Enrich_REACTOME", x = Chemo_Enrich_REACTOME)
        addWorksheet(wb, "Chemo_Enrich_GO")
        writeData(wb, sheet = "Chemo_Enrich_GO", x = Chemo_Enrich_GO)
        saveWorkbook(wb, paste("Chemo_target_enriched",Fold,"Allmol.xlsx", sep = '_'), overwrite = TRUE)


    }


    setwd(Output_file_path)




}
