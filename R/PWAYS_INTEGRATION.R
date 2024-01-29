PWAYS_INTEGRATION = function(Output_file_path,Target_list_path, ORGANISM = 'Mouse', DirPipeline_Data = DirData){

    library(openxlsx)
    library(dplyr)

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

    # CREATE SUMMARY SPREAD SHEETS FOR ENRICHMENT
    FolderType = list.files("../PWAY_ENRICHMENT")

    #KEGG
    for(Fold in FolderType){

        Mols = list.files(paste('../PWAY_ENRICHMENT',Fold, sep = '/'))
        EnrichedPathwaysDf_kegg_All = data.frame("ID" = character(),"Description" = character(), "GeneRatio" = character(),   "BgRatio"=character() , "pvalue" = numeric() ,"p.adjust" = numeric(),  "qvalue" = numeric(), "geneID" = character() , "Count" = integer() ,"Target"= character() ,  "comp" = character() ,"mol" = character())

        for(mol in Mols){

            Filestemp = list.files(paste('../PWAY_ENRICHMENT',Fold,mol, sep = '/'))
            EnrichPwayFile = Filestemp[grep("EnrichPathways_AllGeneSets",Filestemp)]
            load(paste('../PWAY_ENRICHMENT',Fold,mol,EnrichPwayFile, sep = '/'))
            #print(nrow(EnrichedPathwaysDf_kegg))
            EnrichedPathwaysDf_kegg$mol = mol
            EnrichedPathwaysDf_kegg_All = rbind(EnrichedPathwaysDf_kegg_All, EnrichedPathwaysDf_kegg)


        }

        EnrichedPathwaysDf_kegg_All$CompTemp = sapply(EnrichedPathwaysDf_kegg_All$comp, function(x){str_split(x, 'VS')[[1]][1]})
        EnrichedPathwaysDf_kegg_All$NameRef = paste(EnrichedPathwaysDf_kegg_All$comp, EnrichedPathwaysDf_kegg_All$Description, sep = '_')
        EnrichedPathwaysDf_kegg_All$Namesheet = paste(EnrichedPathwaysDf_kegg_All$CompTemp, EnrichedPathwaysDf_kegg_All$Description, sep = '_')
        EnrichedPathwaysDf_kegg_All$Namesheet = sapply(EnrichedPathwaysDf_kegg_All$Namesheet, function(x){substr(x,1,28)})
        seq <- 1:nrow(EnrichedPathwaysDf_kegg_All)
        seq = seq %% 17
        EnrichedPathwaysDf_kegg_All$Namesheet = paste(EnrichedPathwaysDf_kegg_All$Namesheet,seq, sep = '')
        EnrichedPathwaysDf_kegg_All$CompTemp = NULL

        # Create excel workbook
        wb <- createWorkbook()
        addWorksheet(wb, "EnrichedPways_kegg")
        writeData(wb, sheet = "EnrichedPways_kegg", x = EnrichedPathwaysDf_kegg_All)
        # Add the data specific data frames for each pathway
        row = 2
        for(nameid in EnrichedPathwaysDf_kegg_All$NameRef){

            namesheet = EnrichedPathwaysDf_kegg_All$Namesheet[EnrichedPathwaysDf_kegg_All$NameRef == nameid]
            mol = EnrichedPathwaysDf_kegg_All$mol[EnrichedPathwaysDf_kegg_All$NameRef == nameid]
            load(paste('../PWAY_ENRICHMENT',Fold,mol,EnrichPwayFile, sep = '/'))

            if(nameid %in% names(DiffGenesPways_kegg)){
                addWorksheet(wb, namesheet)
                writeData(wb, sheet =  namesheet, x = DiffGenesPways_kegg[[nameid]],rowNames = TRUE)
                writeFormula(wb, "EnrichedPways_kegg",
                             startRow = row, startCol = 14,
                             x = makeHyperlinkString(sheet = namesheet, row = 1, col = 1, text = namesheet)
                )

                writeFormula(wb, namesheet,
                             startRow = 1, startCol = 10,
                             x = makeHyperlinkString(sheet = 'EnrichedPways_kegg', row = 1, col = 1, text = "Back to main")
                )
            }
            row = row+1
        }

        saveWorkbook(wb, paste("EnrichPathways_Kegg",Fold,"_Allmol.xlsx"), overwrite = TRUE)
        save(EnrichedPathwaysDf_kegg_All, file =paste("EnrichPathways_Kegg",Fold,"_Allmol.RData") )









    }

    #REACTOME
    for(Fold in FolderType){

        Mols = list.files(paste('../PWAY_ENRICHMENT',Fold, sep = '/'))
        EnrichedPathwaysDf_Reactome_All = data.frame("ID" = character(),"Description" = character(), "GeneRatio" = character(),   "BgRatio"=character() , "pvalue" = numeric() ,"p.adjust" = numeric(),  "qvalue" = numeric(), "geneID" = character() , "Count" = integer() ,"Target"= character() ,  "comp" = character() ,"mol" = character())

        for(mol in Mols){

            Filestemp = list.files(paste('../PWAY_ENRICHMENT',Fold,mol, sep = '/'))
            EnrichPwayFile = Filestemp[grep("EnrichPathways_AllGeneSets",Filestemp)]
            load(paste('../PWAY_ENRICHMENT',Fold,mol,EnrichPwayFile, sep = '/'))
            #print(nrow(EnrichedPathwaysDf_kegg))
            EnrichedPathwaysDf_Reactome$mol = mol
            EnrichedPathwaysDf_Reactome_All = rbind(EnrichedPathwaysDf_Reactome_All, EnrichedPathwaysDf_Reactome)


        }

        EnrichedPathwaysDf_Reactome_All$CompTemp = sapply(EnrichedPathwaysDf_Reactome_All$comp, function(x){str_split(x, 'VS')[[1]][1]})
        EnrichedPathwaysDf_Reactome_All$NameRef = paste(EnrichedPathwaysDf_Reactome_All$comp, EnrichedPathwaysDf_Reactome_All$Description, sep = '_')
        EnrichedPathwaysDf_Reactome_All$Namesheet = paste(EnrichedPathwaysDf_Reactome_All$CompTemp, EnrichedPathwaysDf_Reactome_All$Description, sep = '_')
        EnrichedPathwaysDf_Reactome_All$Namesheet = sapply(EnrichedPathwaysDf_Reactome_All$Namesheet, function(x){substr(x,1,27)})
        seq <- 1:nrow(EnrichedPathwaysDf_Reactome_All)
        seq = seq %% 100
        EnrichedPathwaysDf_Reactome_All$Namesheet = paste(EnrichedPathwaysDf_Reactome_All$Namesheet,seq, sep = '')
        EnrichedPathwaysDf_Reactome_All$CompTemp = NULL

        # Create excel workbook
        wb <- createWorkbook()
        addWorksheet(wb, "EnrichedPways_Reactome")
        writeData(wb, sheet = "EnrichedPways_Reactome", x = EnrichedPathwaysDf_Reactome_All)
        # Add the data specific data frames for each pathway
        row = 2
        for(nameid in EnrichedPathwaysDf_Reactome_All$NameRef){

            namesheet = EnrichedPathwaysDf_Reactome_All$Namesheet[EnrichedPathwaysDf_Reactome_All$NameRef == nameid]
            mol = EnrichedPathwaysDf_Reactome_All$mol[EnrichedPathwaysDf_Reactome_All$NameRef == nameid]
            load(paste('../PWAY_ENRICHMENT',Fold,mol,EnrichPwayFile, sep = '/'))
            if(nameid %in% names(DiffGenesPways_Reactome)){
                addWorksheet(wb, namesheet)
                writeData(wb, sheet =  namesheet, x = DiffGenesPways_Reactome[[nameid]],rowNames = TRUE)
                writeFormula(wb, "EnrichedPways_Reactome",
                             startRow = row, startCol = 14,
                             x = makeHyperlinkString(sheet = namesheet, row = 1, col = 1, text = namesheet)
                )

                writeFormula(wb, namesheet,
                             startRow = 1, startCol = 10,
                             x = makeHyperlinkString(sheet = 'EnrichedPways_Reactome', row = 1, col = 1, text = "Back to main")
                )
            }
            row = row+1
        }

        saveWorkbook(wb, paste("EnrichPathways_Reactome",Fold,"_Allmol.xlsx"), overwrite = TRUE)
        save(EnrichedPathwaysDf_Reactome_All, file =paste("EnrichPathways_Reactome",Fold,"_Allmol.RData") )



    }


    # SUMMARY OF CHEMOPROTEOMICS

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
        Species = 'mmu'
    }else if(ORGANISM == 'Human'){
        KEGG_GeneSet_path = paste(DirPipeline_Data,"KEGG_All_human_231208.RData", sep = "/")
        REACTOME_GeneSet_path = paste(DirPipeline_Data,"REACTOME_All_human_231209.RData", sep = "/")
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


    # TABLE TO ESTIMATE OVER-REPRESENTATION OF A TARGET WITHIN ENRICHED PATHWAYS

    FolderType = list.files("../PWAY_ENRICHMENT")
    for(Fold in FolderType){

        drugs = list.files(paste("../PWAY_ENRICHMENT",Fold, sep = '/' ))

        # KEGG

        # Create Data frame
        colnames_dyn = c()
        for(treat in drugs){
            colnames_t = c(paste(treat,'Tar_in', sep = '_'), paste(treat,'P_val', sep = '_'))
            colnames_dyn = c(colnames_dyn, colnames_t)

        }

        ColNames = c(c('CompCategory','Chemo_Target', 'Tar_in_KEGG'), colnames_dyn)
        Chemo_Enrich_Kegg = setNames(data.frame(matrix(ncol = length(ColNames), nrow = 0)), ColNames)
        Chemo_Enrich_Kegg[] <- lapply(Chemo_Enrich_Kegg, as.character)

        # Load enriched pathways big summary table
        load(paste("EnrichPathways_Kegg",Fold,"_Allmol.RData") )
        # Add modification to table to make the task easier later
        EnrichedPathwaysDf_kegg_All$CompCategory = sapply(EnrichedPathwaysDf_kegg_All$comp, function(x){paste(str_split(x, '_')[[1]][2:length(str_split(x, '_')[[1]])], collapse = '_')})
        # Load KEGG Gene set
        load(KEGG_GeneSet_path)
        Enriched_WithChemo =  EnrichedPathwaysDf_kegg_All %>% filter(Target != '')
        for(compo in unique(Enriched_WithChemo$CompCategory)){

                Enriched_WithChemo_comp =  Enriched_WithChemo %>% filter(CompCategory == compo)
                Targets_t = sapply(Enriched_WithChemo_comp$Target, function(x){str_split(x, ', ')[[1]] })
                Targets_t = unique(unlist(Targets_t))
                for(ChemoTar in Targets_t){

                    Tar_in_KEGG = paste(as.character(nrow(Chemo_Tar_KEGG[[ChemoTar]])), '/', as.character(nrow(GeneSet_Desc)), sep = '')
                    Tar_in_MOLS = c()
                    Pval_MOLS = c()
                    Joint_vect = c()
                    for(molec in drugs){

                        Enriched_comp_mol =  EnrichedPathwaysDf_kegg_All %>% filter(CompCategory == compo & mol == molec)
                        Enriched_comp_mol_WithTar = Is_Tar_in_Pway(TargetCand = ChemoTar,Target_list, EnrichedPways = Enriched_comp_mol, GeneSet)
                        Tar_in_Mol = nrow(Enriched_comp_mol_WithTar %>% filter(WithTarget != ''))
                        No_Pways_Mol = nrow(Enriched_comp_mol)
                        Tar_in_MOLS[paste(molec,'Tar_in', sep = '_')] = paste(as.character(Tar_in_Mol), "/", as.character(No_Pways_Mol),sep = '')
                        # Hypergeometric distribution: Already sort of a cumulative if use phyper
                        x = Tar_in_Mol # P(X<= x)
                        m = nrow(Chemo_Tar_KEGG[[ChemoTar]]) # Number of pathways containing the target
                        n = nrow(GeneSet_Desc) - m # Number of pathways not containing the target
                        k = No_Pways_Mol # Number of Pathways enriched
                        # x-1 cause I want the right tail from X >= x
                        Pval_MOLS[paste(molec,'P_val', sep = '_')] = sprintf("%e",  1 - phyper(x-1,m,n,k,lower.tail = TRUE ))
                        Joint_vect = c(Joint_vect,c(Tar_in_MOLS[paste(molec,'Tar_in', sep = '_')] ,  Pval_MOLS[paste(molec,'P_val', sep = '_')] ))

                    }

                    Chemo_Enrich_Kegg = Chemo_Enrich_Kegg %>% add_row(data.frame(t(c('CompCategory' = compo,'Chemo_Target' = ChemoTar, 'Tar_in_KEGG' = Tar_in_KEGG,Joint_vect))))
                }
            }

        # Convert Pvalue columns in numeric:
        PvalCol = colnames(Chemo_Enrich_Kegg)[grep('P_val', colnames(Chemo_Enrich_Kegg))]
        for(colname in PvalCol){
            Chemo_Enrich_Kegg[[colname]] = as.numeric(Chemo_Enrich_Kegg[[colname]])
        }

        # REACTOME

        # Create Data frame
        colnames_dyn = c()

        for(treat in drugs){
            colnames_t = c(paste(treat,'Tar_in', sep = '_'), paste(treat,'P_val', sep = '_'))
            colnames_dyn = c(colnames_dyn, colnames_t)

        }

        ColNames = c(c('CompCategory','Chemo_Target', 'Tar_in_REACTOME'), colnames_dyn)
        Chemo_Enrich_React = setNames(data.frame(matrix(ncol = length(ColNames), nrow = 0)), ColNames)
        Chemo_Enrich_React[] <- lapply(Chemo_Enrich_React, as.character)

        # Load enriched pathways big summary table
        load(paste("EnrichPathways_Reactome",Fold,"_Allmol.RData"))
        # Add modification to table to make the task easier later
        EnrichedPathwaysDf_Reactome_All$CompCategory = sapply(EnrichedPathwaysDf_Reactome_All$comp, function(x){paste(str_split(x, '_')[[1]][2:length(str_split(x, '_')[[1]])], collapse = '_')})
        # Load KEGG Gene set
        load(REACTOME_GeneSet_path)
        Enriched_WithChemo_React =  EnrichedPathwaysDf_Reactome_All %>% filter(Target != '')
        for(compo in unique(Enriched_WithChemo_React$CompCategory)){

            Enriched_WithChemo_React_comp =  Enriched_WithChemo_React %>% filter(CompCategory == compo)
            Targets_t = sapply(Enriched_WithChemo_React_comp$Target, function(x){str_split(x, ', ')[[1]] })
            Targets_t = unique(unlist(Targets_t))
            for(ChemoTar in Targets_t){

                Tar_in_React = paste(as.character(nrow(Chemo_Tar_REACTOME[[ChemoTar]])), '/', as.character(nrow(GeneSet_Desc)), sep = '')
                Tar_in_MOLS = c()
                Pval_MOLS = c()
                Joint_vect = c()
                for(molec in drugs){

                    Enriched_comp_mol_React =  EnrichedPathwaysDf_Reactome_All %>% filter(CompCategory == compo & mol == molec)
                    Enriched_comp_mol_WithTar_React = Is_Tar_in_Pway(TargetCand = ChemoTar,Target_list, EnrichedPways = Enriched_comp_mol_React, GeneSet)
                    Tar_in_Mol = nrow(Enriched_comp_mol_WithTar_React %>% filter(WithTarget != ''))
                    No_Pways_Mol = nrow(Enriched_comp_mol_React)
                    Tar_in_MOLS[paste(molec,'Tar_in', sep = '_')] = paste(as.character(Tar_in_Mol), "/", as.character(No_Pways_Mol),sep = '')
                    # Hypergeometric distribution: Already sort of a cumulative if use phyper
                    x = Tar_in_Mol # P(X<= x)
                    m = nrow(Chemo_Tar_REACTOME[[ChemoTar]]) # Number of pathways containing the target
                    n = nrow(GeneSet_Desc) - m # Number of pathways not containing the target
                    k = No_Pways_Mol # Number of Pathways enriched
                    # x-1 cause I want the right tail from X >= x
                    Pval_MOLS[paste(molec,'P_val', sep = '_')] = sprintf("%e",  1 - phyper(x-1,m,n,k,lower.tail = TRUE ))
                    Joint_vect = c(Joint_vect,c(Tar_in_MOLS[paste(molec,'Tar_in', sep = '_')] ,  Pval_MOLS[paste(molec,'P_val', sep = '_')] ))

                }

                Chemo_Enrich_React = Chemo_Enrich_React %>% add_row(data.frame(t(c('CompCategory' = compo,'Chemo_Target' = ChemoTar, 'Tar_in_REACTOME' = Tar_in_React,Joint_vect))))
            }
        }

        # Convert Pvalue columns in numeric:
        PvalCol = colnames(Chemo_Enrich_React)[grep('P_val', colnames(Chemo_Enrich_React))]
        for(colname in PvalCol){
            Chemo_Enrich_React[[colname]] = as.numeric(Chemo_Enrich_React[[colname]])
        }

        # Create excel workbook
        wb <- createWorkbook()
        addWorksheet(wb, "Chemo_Enrich_Kegg")
        writeData(wb, sheet = "Chemo_Enrich_Kegg", x = Chemo_Enrich_Kegg)
        addWorksheet(wb, "Chemo_Enrich_React")
        writeData(wb, sheet = "Chemo_Enrich_React", x = Chemo_Enrich_React)
        saveWorkbook(wb, paste("Chemo_target_enriched",Fold,"_Allmol.xlsx"), overwrite = TRUE)


    }







}
