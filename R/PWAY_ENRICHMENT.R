# DEG per time point and Pathway Enrichment
PWAY_ENRICHMENT = function(Output_file_path, WITH_REVERSE = TRUE, REVERSE_DRUG = 'Sham', REVERSE_TIME = "72", padjval = 0.2, LogFoldThrs = 1,Pway_qvalThrs = 0.2, ORGANISM = 'Mouse', Target_list_path, DirPipeline_Data, PARAllEL = TRUE, Partition = 'compute' , GeneDescription_path){
    library(ggvenn)
    library(ggplot2)
    library(ggrepel)
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    library(clusterProfiler)
    library(pathview)
    library(DESeq2)
    library(eulerr)
    library(org.Hs.eg.db)
    library(stringr)
    library(clustermq)
    library(foreach)
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(rstatix)
    library(ggpubr)



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

    # Function to add description and gene plots of regulated genes in a pathway
    Gene_Expr_plot_and_Descrip = function(ti,NormMatrix,Metadata, Gene_comb, ControlName, Df_mol_filt, comp, GeneDescription, mol){

        gids = unlist(str_split(Gene_comb[ti], '/'))
        GeneList = Df_mol_filt %>% filter(entrez %in% gids)
        GeneList = GeneList$genes
        NormMatrix_filt = as.data.frame(NormMatrix[GeneList,])
        NormMatrix_filt$Gene = rownames(NormMatrix_filt)
        NormMatrix_Long = NormMatrix_filt %>% pivot_longer(!Gene, names_to = "Sample_name", values_to = "Norm_Expr") %>% as.data.frame()
        NormMatrix_Long$Gene = factor(NormMatrix_Long$Gene)
        NormMatrix_Long$time = Metadata$Rna_collection_time_hrs[match(NormMatrix_Long$Sample_name, Metadata$Sample_name)]
        NormMatrix_Long$time = factor(NormMatrix_Long$time, levels = sort(unique(NormMatrix_Long$time)))
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

        png(filename=paste('Genes', comp, gsub("/", "_", Gene_comb[ti]), '.png',sep = '_'), width = 1200, height = 800, res= 100)
        print(bxp)
        dev.off()

        # Add gene data frame with description

        Gene_info_temp = Df_mol_filt %>% filter(entrez %in% gids) %>% dplyr::select("log2FoldChange", "pvalue", "padj", "genes", "entrez")

        # Gene description data frame
        Description = c()
        for(gene in Gene_info_temp$genes){

            for(i in 1:nrow(GeneDescription)){

                if(gene %in% GeneDescription$Gene.Names_mod[i][[1]]){
                    Description[gene] = GeneDescription$LongDesc[i]

                    break
                }

            }
        }

        Gene_info_temp$Description = Description[match(Gene_info_temp$genes, names(Description))]
        DiffGenesPwaysDf = list()
        DiffGenesPwaysDf[[paste('Genes', comp, gsub("/", "_", Gene_comb[ti]), '.png',sep = '_')]] = Gene_info_temp
        DiffGenesPwaysDf[['Names']] = paste('Genes', comp, gsub("/", "_", Gene_comb[ti]), '.png',sep = '_')

        return(DiffGenesPwaysDf)



    }

    # Function to plot Kegg diagrams
    Plot_kegg = function(ti,enrichedPathways, Df_mol_filt, comp, Problematic_Pways) {

        curKegg = enrichedPathways$ID[ti]
        curKegg_Desc = enrichedPathways$Description[ti]
        gids = unlist(str_split(enrichedPathways$geneID[ti], '/'))
        GeneFold = Df_mol_filt$log2FoldChange[which(Df_mol_filt$entrez %in% gids)]
        names(GeneFold) = Df_mol_filt$entrez[which(Df_mol_filt$entrez %in% gids)]
        DiffGenesPways_kegg_t= Df_mol_filt %>% filter(entrez %in% gids) %>% dplyr::select("log2FoldChange", "pvalue", "padj", "genes", "entrez")
        NamePathway = paste(comp,curKegg_Desc,sep='_' )
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

        DiffGenesPways_kegg_list = list("DiffGenesPways_kegg_t" = DiffGenesPways_kegg_t,"NamePathway" = NamePathway )
        return(DiffGenesPways_kegg_list)
    }
    # Cluster temp file
    ClusterTemp = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/slurmMqBoris.tmpl"


    # Checking parameters if FILTER_REVERSE = FALSE does not matter what REVERSE_CONTROL is, cause it is not used. If FILTER_REVERSE = TRUE, the length of REVERSE_CONTROL
    # should be bigger than zero.
    # Data bases for pathway enrichment
    Gene_set_names = c('KEGG', 'REACTOME', 'GO')
    if(ORGANISM == 'Mouse'){
        # Mouse
        KEGG_GeneSet_path = paste(DirPipeline_Data,"KEGG_All_mouse_230922.RData", sep = "/")
        REACTOME_GeneSet_path = paste(DirPipeline_Data,"REACTOME_All_mouse_230924.RData", sep = "/")
        GO_GeneSet_path = paste(DirPipeline_Data,"GO_All_mouse_240212.RData", sep = "/")

        Species = 'mmu'
        # Read gene description
        GeneDescription <- read_excel(GeneDescription_path)
        # Re do this step because when we save it as .xlsx, the vector properties of the column Gene.Names_mod are lost
        GeneDescription$Gene.Names_mod = sapply(GeneDescription$Gene.Names,function(x) {str_split(x,' ')[[1]]})
        GeneDescription = GeneDescription %>% filter(!is.na(Function..CC.))

    }else if(ORGANISM == 'Human'){
        KEGG_GeneSet_path = paste(DirPipeline_Data,"KEGG_All_human_231208.RData", sep = "/")
        REACTOME_GeneSet_path = paste(DirPipeline_Data,"REACTOME_All_human_231209.RData", sep = "/")
        GO_GeneSet_path = paste(DirPipeline_Data,"GO_All_human_240212.RData", sep = "/")
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

    # Download normalized matrix
    Data_Norm = "./DESEQ_NORM_QCNORM/DESeq_Norm.RData"
    load(Data_Norm)


    # Create Folder
    setwd(Output_file_path)
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

                # Define Universe of genes for future enrichment

                Univ_genes = rownames(DiffGenes[[1]])
                Univ_genes_entrez <- mapIds(org.Mm.eg.db, keys=Univ_genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
                naIdx = which(is.na(Univ_genes_entrez))
                if (length(naIdx) > 0) {
                    Univ_genes_entrez =  Univ_genes_entrez[-naIdx]
                }


                # List of comparisons or pathway enrichment
                Listdrug = names(DiffGenes)[grep(drug,names(DiffGenes))]

                DiffGenesPways_KEGG = list()
                DiffGenesPways_REACTOME = list()
                DiffGenesPways_GO = list()

                DiffUniqueGenesPwaysDf_Mol_KEGG = list()
                DiffUniqueGenesPwaysDf_Mol_REACTOME = list()
                DiffUniqueGenesPwaysDf_Mol_GO = list()

                EnrichedPathwaysDf_KEGG = data.frame("ID" = character(), "Description" = character(), "GeneRatio" = character, "BgRatio" = character,
                                                     "pvalue" = numeric(), "p.adjust" = numeric(), "qvalue" = numeric(), "geneID" = character(),
                                                     "Count" = integer(), "Target" = character(), "comp" = character() )

                EnrichedPathwaysDf_REACTOME = data.frame("ID" = character(), "Description" = character(), "GeneRatio" = character, "BgRatio" = character,
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


                    for(Gene_set_name in Gene_set_names){

                        # Create Folder for gene set
                        Name_folder_L3 = Gene_set_name
                        dir.create(Name_folder_L3)
                        setwd(Name_folder_L3)


                        print(paste('Enrich', Gene_set_name, sep = '_'))

                        if(Gene_set_name == 'KEGG'){

                            load(KEGG_GeneSet_path)
                            try(res <- enricher(Df_mol_filt$entrez, TERM2GENE = GeneSet, TERM2NAME=GeneSet_Desc, universe = Univ_genes_entrez))
                        }else if(Gene_set_name == 'REACTOME'){

                            load(REACTOME_GeneSet_path)
                            try(res <- enricher(Df_mol_filt$entrez, TERM2GENE = GeneSet, TERM2NAME=GeneSet_Desc, universe = Univ_genes_entrez))

                        }else if(Gene_set_name == 'GO'){


                            load(GO_GeneSet_path)
                            try(res <- enricher(Df_mol_filt$entrez, TERM2GENE = GeneSet, TERM2NAME=GeneSet_Desc, universe = Univ_genes_entrez))

                        }



                        if(!is.null(res)){

                            enrichedPathways = res@result %>% filter(p.adjust < Pway_qvalThrs & Count > 1)
                            if(nrow(enrichedPathways) > 0){
                            #if(nrow(summary(res)) >1 & max(summary(res)$Count) > 1){

                                # Fixing notation
                                resultTemp = res@result
                                if(ORGANISM == 'Mouse' & Gene_set_name == 'KEGG'){

                                   resultTemp$Description = sapply(resultTemp$Description, function(x){str_split(x,' - Mus')[[1]][1]})
                                   res@result = resultTemp
                                   enrichedPathways$Description = sapply(enrichedPathways$Description, function(x){str_split(x,' - Mus')[[1]][1]})

                                }



                                # Add chemoproteomics candidates to the pathways enriched in the data frame
                                enrichedPathways = Add_Chemo_Enrich(drug, Target_list, enrichedPathways)
                                # Save data frame of enriched pathways
                                enrichedPathways$comp = comp
                                #EnrichedPathwaysDf_KEGG = rbind(EnrichedPathwaysDf_KEGG,enrichedPathways )

                                # Dotplot
                                # png(file="DotPlot_Enrich.png",width = 800, height = 800, res = 100)
                                # print(try(dotplot(res, showCategory=20, font.size = 12) + ggtitle(paste(drug, comp, sep = ' '))))
                                # dev.off()

                                # cnetplot
                                if(ORGANISM == 'Mouse'){
                                    #This data frame for the cnetplot
                                    resx <- setReadable(res, "org.Mm.eg.db", 'ENTREZID')

                                }

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


                                if(PARAllEL == TRUE){

                                    print("Running pathway enrichment in parallel")

                                    TIMEOUT = 1000
                                    NETWORKS = min(nrow(enrichedPathways), 100)

                                    options(
                                        clustermq.scheduler = "slurm",
                                        clustermq.template = ClusterTemp,
                                        clustermq.data.warning=5000 #megabytes
                                    )
                                    register_dopar_cmq(n_jobs=NETWORKS,
                                                       fail_on_error=FALSE,
                                                       verbose=TRUE,
                                                       log_worker=TRUE,
                                                       timeout = TIMEOUT,
                                                       template=list(
                                                           timeout=TIMEOUT, #how long to wait on SLURM side
                                                           memory=5000,
                                                           cores=1,#how many cores to use (to throttle down memory usage),
                                                           partition = 'compute',
                                                           r_path = file.path(R.home("bin"), "R")))


                                    if(Gene_set_name == 'KEGG'){

                                            EnrichRes = foreach (ti = 1:nrow(enrichedPathways),.export=c( 'Plot_kegg'),
                                                                 .packages=c('stringr', 'pathview')) %dopar% {


                                                                     res <- Plot_kegg(ti,enrichedPathways, Df_mol_filt, comp, Problematic_Pways)
                                                                     return(res)

                                                                 }
                                            DiffGenesPways_kegg_temp = lapply(EnrichRes, function(x){x$DiffGenesPways_kegg_t})
                                            names(DiffGenesPways_kegg_temp) = lapply(EnrichRes, function(x){x$NamePathway})

                                            Files_Current = list.files(".")
                                            Log_files = Files_Current[startsWith(Files_Current, "cmq")]
                                            Name_folder_LOGS = "Logs_KEGG"
                                            dir.create(Name_folder_LOGS)
                                            file.copy(file.path(".",Log_files),"./Logs_KEGG")
                                            file.remove(file.path(".",Log_files))

                                            DiffGenesPways_KEGG = c(DiffGenesPways_KEGG,DiffGenesPways_kegg_temp)
                                    }

                                    # Add gene description to data fame and plot gene expression

                                    Unique_gene_comb = unique(enrichedPathways$geneID)

                                    # THIS SECTION NEEDS TO BE MODIFIED WHEN WE HAVE DOSE ALSO VARYING

                                    Control = str_split(str_split(comp,  '_VS_')[[1]][2],'_')[[1]][1]
                                    #Time_coll = str_split(str_split(comp,  '_VS_')[[1]][2],'_')[[1]][2]
                                    Metadata_comp_drug = Metadata %>% filter(Treatment %in% c(drug, Control))
                                    Norm_comp_drug = as.data.frame(NormCounts[,match(Metadata_comp_drug$Sample_name, colnames(NormCounts))])

                                    DiffUniqueGenesPwaysDf_temp = foreach (ti = 1:length(Unique_gene_comb),.export=c( 'Genes_df_desc'),
                                                         .packages=c('stringr', 'tidyr','rstatix', 'ggpubr')) %dopar% {
                                                        res <- Gene_Expr_plot_and_Descrip(ti,NormMatrix = Norm_comp_drug, Metadata = Metadata_comp_drug, Gene_comb = Unique_gene_comb, ControlName = Control,  Df_mol_filt = Df_mol_filt, comp = comp, GeneDescription = GeneDescription, mol = drug)
                                                        return(res)

                                                         }

                                    DiffUniqueGenesPwaysDf = lapply(DiffUniqueGenesPwaysDf_temp, function(x){x[[1]]})
                                    names(DiffUniqueGenesPwaysDf) = unlist(lapply(DiffUniqueGenesPwaysDf_temp, function(x){x[[2]]}))

                                    # Move Log files

                                    Files_Current = list.files(".")
                                    Log_files = Files_Current[startsWith(Files_Current, "cmq")]
                                    Name_folder_LOGS = "Logs_Expr"
                                    dir.create(Name_folder_LOGS)
                                    file.copy(file.path(".",Log_files),"./Logs_Expr")
                                    file.remove(file.path(".",Log_files))


                                }else{

                                    EnrichRes = list()
                                    DiffUniqueGenesPwaysDf_temp = list()

                                    if(Gene_set_name == 'KEGG'){

                                        for(ti in 1:nrow(enrichedPathways)){

                                            print("Running pathway enrichment in serial model")

                                            EnrichRes[[ti]] = Plot_kegg(ti,enrichedPathways, Df_mol_filt, comp, Problematic_Pways)


                                        }

                                        DiffGenesPways_kegg_temp = lapply(EnrichRes, function(x){x$DiffGenesPways_kegg_t})
                                        names(DiffGenesPways_kegg_temp) = lapply(EnrichRes, function(x){x$NamePathway})
                                        DiffGenesPways_KEGG = c(DiffGenesPways_KEGG,DiffGenesPways_kegg_temp)

                                    }

                                    # Add gene description to data fame and plot gene expression

                                    Unique_gene_comb = unique(enrichedPathways$geneID)

                                    # THIS SECTION NEEDS TO BE MODIFIED WHEN WE HAVE DOSE ALSO VARYING

                                    Control = str_split(str_split(comp,  '_VS_')[[1]][2],'_')[[1]][1]
                                    Time_coll = str_split(str_split(comp,  '_VS_')[[1]][2],'_')[[1]][2]
                                    Metadata_comp_drug = Metadata %>% filter(Treatment %in% c(drug, Control))
                                    Norm_comp_drug = as.data.frame(NormCounts[,match(Metadata_comp_drug$Sample_name, colnames(NormCounts))])



                                    for(ti in 1:length(Unique_gene_comb)){

                                        DiffUniqueGenesPwaysDf_temp[[ti]] = Gene_Expr_plot_and_Descrip(ti,NormMatrix = Norm_comp_drug, Metadata = Metadata_comp_drug, Gene_comb = Unique_gene_comb, ControlName = Control,  Df_mol_filt = Df_mol_filt, comp = comp, GeneDescription = GeneDescription, mol = drug)
                                    }

                                    DiffUniqueGenesPwaysDf = lapply(DiffUniqueGenesPwaysDf_temp, function(x){x[[1]]})
                                    names(DiffUniqueGenesPwaysDf) = unlist(lapply(DiffUniqueGenesPwaysDf_temp, function(x){x[[2]]}))


                                }

                                if(Gene_set_name == 'KEGG'){

                                    EnrichedPathwaysDf_KEGG = rbind(EnrichedPathwaysDf_KEGG,enrichedPathways )
                                    DiffUniqueGenesPwaysDf_Mol_KEGG = c(DiffUniqueGenesPwaysDf_Mol_KEGG, DiffUniqueGenesPwaysDf)

                                }else if(Gene_set_name == 'REACTOME'){

                                    EnrichedPathwaysDf_REACTOME = rbind(EnrichedPathwaysDf_REACTOME,enrichedPathways )
                                    DiffUniqueGenesPwaysDf_Mol_REACTOME = c(DiffUniqueGenesPwaysDf_Mol_REACTOME, DiffUniqueGenesPwaysDf)

                                }else if(Gene_set_name == 'GO'){

                                    EnrichedPathwaysDf_GO = rbind(EnrichedPathwaysDf_GO,enrichedPathways )
                                    DiffUniqueGenesPwaysDf_Mol_GO = c(DiffUniqueGenesPwaysDf_Mol_GO, DiffUniqueGenesPwaysDf)
                                }



                            }else{
                                print(paste(Gene_set_name, "enrichment is not very reliable, therefore not saved", sep = ' '))
                            }
                        }

                        setwd("..")
                    }

                    # #REACTOME
                    #
                    # Name_folder_L3 = "REACTOME"
                    # dir.create(Name_folder_L3)
                    # setwd(Name_folder_L3)
                    #
                    #
                    # print('Enrich REACTOME')
                    # # Get Reactome from mouse
                    #
                    # load(REACTOME_GeneSet_path)
                    # try(res <- enricher(Df_mol_filt$entrez, TERM2GENE = GeneSet,TERM2NAME = GeneSet_Desc, universe = Univ_genes_entrez ))
                    #
                    # if(!is.null(res)){
                    #
                    #   enrichedPathways = res@result %>% filter(p.adjust < Pway_qvalThrs & Count > 1)
                    #   if(nrow(enrichedPathways) > 0){
                    #     #if(nrow(summary(res)) >1 & max(summary(res)$Count) > 1){
                    #     #enrichedPathways = res@result %>% filter(qvalue < Pway_qvalThrs)
                    #
                    #     # Add chemoproteomics candidates to the pathways enriched in the data frame
                    #     enrichedPathways = Add_Chemo_Enrich(drug, Target_list, enrichedPathways)
                    #
                    #     # Save enriched pathways
                    #     enrichedPathways$comp = comp
                    #     EnrichedPathwaysDf_Reactome = rbind(EnrichedPathwaysDf_Reactome,enrichedPathways )
                    #     #save(enrichedPathways, file = paste("EnrichPathways","Pway_qvalThrs", Pway_qvalThrs,".RData", sep = "_"))
                    #
                    #
                    #
                    #     # Plots
                    #     Dot = dotplot(res, showCategory=20, font.size = 12) + ggtitle(paste(drug, comp, sep = ' '))
                    #     # Dotplot
                    #     png(file="DotPlot_Enrich.png",width = 800, height = 800, res = 100)
                    #     print(Dot)
                    #     dev.off()
                    #
                    #     # cnetplot
                    #     # Adding symbol
                    #     if(ORGANISM == 'Mouse'){
                    #         # This data frame for the cnetplot
                    #         resx <- setReadable(res, "org.Mm.eg.db", 'ENTREZID')
                    #
                    #     }
                    #     genelist = Df_mol_filt$log2FoldChange
                    #     names(genelist) = Df_mol_filt$entrez
                    #
                    #     png(file="CnetPlot_Enrich.png",width = 800, height = 800, res = 100)
                    #     if(ORGANISM == "Mouse"){
                    #         print(try(cnetplot(resx, foldChange=genelist, showCategory = 10)))
                    #     }else{
                    #
                    #         print(try(cnetplot(res, foldChange=genelist, showCategory = 10)))
                    #     }
                    #     dev.off()
                    #
                    #
                    #     Problematic_Pways = c()
                    #     for (ti in 1:length(enrichedPathways$ID)){
                    #         if (!(enrichedPathways$ID[ti] %in% Problematic_Pways)){
                    #             curKegg = enrichedPathways$ID[ti]
                    #             print(curKegg)
                    #             curKegg_Desc = enrichedPathways$Description[ti]
                    #             gids = unlist(str_split(enrichedPathways$geneID[ti], '/'))
                    #             Gene_info_temp = Df_mol_filt %>% filter(entrez %in% gids) %>% dplyr::select("log2FoldChange", "pvalue", "padj", "genes", "entrez")
                    #
                    #             # Gene description data frame
                    #             Description = c()
                    #             for(gene in Gene_info_temp$genes){
                    #
                    #                 for(i in 1:nrow(GeneDescription)){
                    #
                    #                     if(gene %in% GeneDescription$Gene.Names_mod[i][[1]]){
                    #                         Description[gene] = GeneDescription$LongDesc[i]
                    #
                    #                         break
                    #                     }
                    #
                    #                 }
                    #             }
                    #
                    #             Gene_info_temp$Description = Description[match(Gene_info_temp$genes, names(Description))]
                    #             DiffGenesPways_Reactome[[paste(drug, comp,curKegg_Desc,sep='_' )]] = Gene_info_temp
                    #
                    #         }
                    #     }
                    #
                    #   }else{
                    #
                    #       print("REACTOME enrichment is not very reliable, therefore not saved")
                    #
                    #    }
                    # }
                    #
                    # setwd("..")



                    # GO Enrichment

                    # Name_folder_L3 = "GO"
                    # dir.create(Name_folder_L3)
                    # setwd(Name_folder_L3)
                    #
                    # print('Enrich GO')
                    #
                    # if(ORGANISM == 'Mouse'){
                    #
                    #     try(res <- enrichGO(Df_mol_filt$entrez, OrgDb = "org.Mm.eg.db", ont="ALL",readable=TRUE, qvalueCutoff = Pway_qvalThrs, universe = Univ_genes_entrez))
                    #
                    # }else if(ORGANISM == 'Human'){
                    #
                    #
                    #     try(res <- enrichGO(Df_mol_filt$entrez, OrgDb = 'org.Hs.eg.db', ont="ALL",readable=TRUE, qvalueCutoff = Pway_qvalThrs, universe = Univ_genes_entrez))
                    # }else{
                    #
                    #     stop("Select a valid organism")
                    # }
                    #
                    # if(!is.null(res)){
                    #
                    #    enrichedPathways = res@result %>% filter(p.adjust < Pway_qvalThrs & Count > 1)
                    #    if(nrow(enrichedPathways) > 0){
                    #
                    #     # if(nrow(summary(res)) >1 & max(summary(res)$Count) > 1){
                    #     #
                    #     # enrichedPathways = res@result %>% filter(p.adjust < Pway_qvalThrs)
                    #
                    #     # Add chemoproteomics candidates to the pathways enriched in the data frame
                    #     enrichedPathways = Add_Chemo_Enrich(drug, Target_list, enrichedPathways)
                    #
                    #     # Save enriched pathways
                    #     enrichedPathways$comp = comp
                    #     EnrichedPathwaysDf_GO = rbind(EnrichedPathwaysDf_GO,enrichedPathways )
                    #     #save(enrichedPathways, file = paste("EnrichPathways","Pway_qvalThrs", Pway_qvalThrs,".RData", sep = "_"))
                    #
                    #
                    #     # Plots
                    #
                    #     # Dotplot
                    #     png(file="DotPlot.png", width = 800, height = 800, res = 100)
                    #     try(print(dotplot(res, showCategory=20, font.size = 10, split = "ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")))
                    #     dev.off()
                    #
                    #     # cnetplot
                    #     # Adding symbol
                    #     if(ORGANISM == 'Mouse'){
                    #         # This data frame for the cnetplot
                    #         resx <- setReadable(res, "org.Mm.eg.db", 'ENTREZID')
                    #
                    #     }
                    #     genelist = Df_mol_filt$log2FoldChange
                    #     names(genelist) = Df_mol_filt$entrez
                    #
                    #     png(file="CnetPlot_Enrich.png",width = 800, height = 800, res = 100)
                    #     if(ORGANISM == "Mouse"){
                    #         print(try(cnetplot(resx, foldChange=genelist, showCategory = 10)))
                    #     }else{
                    #
                    #         print(try(cnetplot(res, foldChange=genelist, showCategory = 10)))
                    #     }
                    #     dev.off()
                    #
                    #     Problematic_Pways = c()
                    #     for (ti in 1:length(enrichedPathways$ID)){
                    #         if (!(enrichedPathways$ID[ti] %in% Problematic_Pways)){
                    #             curKegg = enrichedPathways$ID[ti]
                    #             print(curKegg)
                    #             curKegg_Desc = enrichedPathways$Description[ti]
                    #             gids = unlist(str_split(enrichedPathways$geneID[ti], '/'))
                    #             GeneFold = Df_mol_filt$log2FoldChange[which(Df_mol_filt$Gene %in% gids)]
                    #             names(GeneFold) = Df_mol_filt$entrez[which(Df_mol_filt$Gene %in% gids)]
                    #             print(GeneFold)
                    #             DiffGenesPways_GO[[paste(comp,curKegg_Desc,sep='_' )]] = Df_mol_filt %>% filter(genes %in% gids) %>% dplyr::select("log2FoldChange", "pvalue", "padj", "genes", "entrez")
                    #
                    #         }
                    #     }
                    #    }else{
                    #
                    #        print("GO enrichment is not very reliable, therefore not saved")
                    #
                    #   }
                    # }

                    setwd("..")

                }


                save(EnrichedPathwaysDf_KEGG,EnrichedPathwaysDf_REACTOME,EnrichedPathwaysDf_GO, DiffGenesPways_KEGG, DiffGenesPways_REACTOME, DiffGenesPways_GO ,
                     DiffUniqueGenesPwaysDf_Mol_KEGG,DiffUniqueGenesPwaysDf_Mol_REACTOME,DiffUniqueGenesPwaysDf_Mol_GO,
                     file = paste('EnrichPathways_AllGeneSets', 'FILTER_REVERSE', FILTER_REVERSE,'.RData', sep = '_'))
                setwd('..')

            }

            setwd('..')
    }

    setwd(Output_file_path)
}








