# DESeq normalization and DEG
DESEQ_NORM = function(Output_file_path, QCNORM = "PRE_QCNORM", outliers_path, CoarseConditions, METHOD_NORM = 'Standard', VST_FILTER = "VST_ON",  SVD_FILTER = 'SVD_OFF', PlotPCA = 'PCA_PLOT', Control_Neg_PCA = "", Control_Pos_PCA = "", MEM_MB = 1600, GENE_ENSEM = TRUE, ORGANISM = 'Human', Batch_variable = 'Plate.id' ){

    library(DESeq2)
    #library(readxl)
    library(dplyr)
    library(BiocParallel)
    library(parallel)
    library(PCAtools)
    library(foreach)
    library(doParallel)
    library(openxlsx)
    nCores = detectCores()
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    library(clustermq)
    library(org.Hs.eg.db)
    library("cowplot")
    library('pheatmap')
    library("GGally")
    library("sva")



    source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/Functions_Invivo.R")
    clusterTemp = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/slurmMqBoris.tmpl"

    # Import data
    setwd(Output_file_path)
    Data_file = "./PRE_FILTERING/Prefilter_Data.xlsx"
    Count = read.xlsx(xlsxFile = Data_file, sheet = "Count", rowNames= TRUE)
    Metadata = read.xlsx(xlsxFile = Data_file, sheet = "Metadata")

    # Read customize outliers
    if(file.exists(outliers_path)){

        #outliers_cust <- read.csv(outliers_path, header = TRUE)
        outliers_cust = read.xlsx(outliers_path, colNames = TRUE)
        outliers_cust = outliers_cust[[1]]
    }else{

        print("Path to costumized outliers does not exist or is null")
        outliers_cust = c()
    }

    # Add outlier feature to metadata
    Metadata$Outliers = 'No'
    Metadata$OutliersNames = ''
    if(length(outliers_cust)>0){

        Metadata$Outliers[Metadata$Sample_name %in% outliers_cust] = 'Yes'
        Metadata$OutliersNames = ''
        Metadata$OutliersNames[Metadata$Outliers=='Yes'] = Metadata$Sample_name[Metadata$Outliers=='Yes']

    }



    # Evaluate condition depending on whether the function is call before or after QC_POSTNORNMALIZATION
    if(QCNORM == "POST_QCNORM"){
        outliers = outliers_cust


        # Create Folder

        Name_folder =  "DESEQ_NORM_QCNORM"
        dir.create(Name_folder)
        setwd(Name_folder)

        # Filter for outliers

        if(length(outliers)>0){
            print('Getting rid of outliers')

            Metadata = Metadata %>% filter(!Sample_name %in% outliers)
            Count = Count[ ,match(Metadata$Sample_name,colnames(Count) ) ]
            print(dim(Count))
        }else{
            print('No outliers list provided')
        }



    }else if(QCNORM == "PRE_QCNORM"){

        Name_folder =  "DESEQ_NORM"
        dir.create(Name_folder)
        setwd(Name_folder)

        outliers = outliers_cust

    }else{

        stop("Select a valid QC_NORM option")
    }



    # DESeq normalization
    #Create coarse Condition feature and useful Metadata
    Metadata$CoarseCondition = Metadata[[CoarseConditions[1]]]

    for(j in 2:length(CoarseConditions)){
        Metadata$CoarseCondition = paste(Metadata[["CoarseCondition"]],Metadata[[CoarseConditions[j]]], sep = "_")
    }

    # Convert from ensemble to symbol
    if(GENE_ENSEM == TRUE){
        if(ORGANISM == 'Human'){
        # Convert from Ensemble to symbol
        geneSymbol <- mapIds(org.Hs.eg.db,keys=rownames(Count),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
        Count$genesymbol = geneSymbol[match(rownames(Count), names(geneSymbol))]
        Count = Count[!is.na(Count$genesymbol),]

        Count <- aggregate(. ~ genesymbol, data = Count, FUN = sum)
        rownames(Count) = Count$genesymbol
        Count$genesymbol = NULL
      }
    }


    # Normalization
    if (METHOD_NORM == 'Standard'){
        print('Running DESEQ Standard normalization...')
        ##create a DESeq object

        Design = "~ CoarseCondition"
        # To be modified if more than one factor is desired to be adjusted for
        if(length(Batch_variable) == 1 ){

            Design = paste(Design, Batch_variable, sep = ' + ' )



        }

        dds <- DESeqDataSetFromMatrix(countData = Count, colData = Metadata, design = as.formula(Design))
        startTime <- Sys.time()
        deseqObj = DESeq(dds,
                         parallel = TRUE,
                         fitType = "parametric",
                         BPPARAM=MulticoreParam(detectCores() - 2))
        if(VST_FILTER == "VST_ON"){
            print("Using Vst ")
            NormCounts = getVarianceStabilizedData(deseqObj)
            endTime <- Sys.time()
            print(endTime - startTime)
        }else if(VST_FILTER == "VST_OFF"){
            print("Not using vst")
            NormCounts = counts(deseqObj, normalized=TRUE)

        }else{

            stop('Choose VST_FILTER method valid')
        }

    }else if(METHOD_NORM == 'Standard_Parallel'){

        print('Running DESEQ Standard normalization in parallel...')
        NJOBS = 100
        TIMEOUT = 10000
        MEMORY = MEM_MB

        options(
            clustermq.scheduler = "slurm",
            clustermq.template = clusterTemp,
            clustermq.data.warning=5000 #megabytes
        )
        register(DoparParam())
        register_dopar_cmq(n_jobs=NJOBS, memory=MEMORY, pkgs="BiocParallel", export=list(
            .bpworker_EXEC=BiocParallel:::.bpworker_EXEC,
            .log_buffer_get=BiocParallel:::.log_buffer_get,
            #.log_data=BiocParallel:::.log_data,
            .log_buffer_init=BiocParallel:::.log_buffer_init,
            .VALUE=BiocParallel:::.VALUE
        ),
        template=list(
            timeout=TIMEOUT, #how long to wait on SLURM side
            memory=5000,
            cores=1,#how many cores to use (to throttle down memory usage),
            partition = 'himem',
            r_path = file.path(R.home("bin"), "R")))


            Design = "~ CoarseCondition"
            # To be modified if more than one factor is desired to be adjusted for
            if(length(Batch_variable) == 1 ){

                Design = paste(Design, Batch_variable, sep = ' + ' )



            }


            dds <- DESeqDataSetFromMatrix(countData = Count, colData = Metadata, design = as.formula(Design))
            startTime <- Sys.time()
            deseqObj = DESeq(dds,
                             parallel = TRUE,
                             fitType = "parametric",
                             BPPARAM=bpparam())
            if(VST_FILTER == "VST_ON"){
                print("Using Vst ")
                NormCounts = getVarianceStabilizedData(deseqObj)
                endTime <- Sys.time()
                print(endTime - startTime)
            }else if(VST_FILTER == "VST_OFF"){
                print("Not using vst")
                NormCounts = counts(deseqObj, normalized=TRUE)

            }else{

                stop('Choose VST_FILTER method valid')
            }



    }else {

        stop('Choose a valid method for normalization')
    }

   # SVD filter
    if(SVD_FILTER == 'SVD_ON'){

        print('Truncated SVD ON')
        # Get significant PCA
        hornRes = parallelPCA(
            NormCounts,
            max.rank = min(nrow(NormCounts), ncol(NormCounts)),
            niters = 200,
            center = TRUE,
            scale=TRUE,
            threshold = 0.1,
            transposed = FALSE,
            BPPARAM = DoparParam() #not parallelizing since we are now parallelizing top level
        )

        subDf =t(NormCounts)
        NormCounts_scaled = apply(subDf, 2, scale)
        rownames(NormCounts_scaled) = rownames(subDf)
        S = svd(NormCounts_scaled)
        #M  = S$u %*% diag(S$d) %*% t(S$v)
        #Max = max(abs(M - NormCounts_scaled))

        # Let's take only significant components of the svd
        Sig_d = S$d
        Sig_d[(hornRes$n+1) : length(Sig_d)] = 0
        NormCounts_scaled_trunc =  S$u %*% diag(Sig_d) %*% t(S$v)
        rownames(NormCounts_scaled_trunc) = rownames(NormCounts_scaled)
        colnames(NormCounts_scaled_trunc) = colnames(NormCounts_scaled)
        # Just checking
        #MT =  S$u[,1:hornRes$n] %*% diag(S$d[1:hornRes$n]) %*% t(S$v[,1:hornRes$n])
        #Max = max(abs(MT - NormCounts_scaled_trunc))
        NormCounts = t(NormCounts_scaled_trunc)

    }

    # Correct for batch
    if(Batch_variable != ''){
        pheno = Metadata %>% dplyr::select(Sample_name, Plate.id, CoarseCondition) %>% column_to_rownames(var = "Sample_name")
        pheno[[Batch_variable]] = as.factor(pheno[[Batch_variable]])
        modcombat = model.matrix(~CoarseCondition, data = pheno)
        NormCounts_corr = ComBat(dat= NormCounts, batch= pheno[[Batch_variable]], mod=modcombat, par.prior=TRUE)
        NormCounts = NormCounts_corr
    }


    #Save RData object
    save(deseqObj, NormCounts,Metadata, file = 'DESeq_Norm.RData')
    # Save excel file
    wb <- createWorkbook()
    addWorksheet(wb, "NormCounts")
    addWorksheet(wb, "Metadata")
    writeData(wb, sheet = "NormCounts", x = NormCounts, rowNames = TRUE)
    writeData(wb, sheet = "Metadata", x = Metadata)
    saveWorkbook(wb, "Norm_Data.xlsx", overwrite = TRUE)


    if(PlotPCA == 'PCA_PLOT'){


        Name_folder_PCA =  "PCA_PLOTS"
        dir.create(Name_folder_PCA)
        setwd(Name_folder_PCA )

        # Create indication On feature in Metadata from Indication_induction_time(hrs) or from "Stimulant_used.
        # Only one of those should be present depending on whether cell line or in vivo metadata

        if("Indication_induction_time_hrs" %in% colnames(Metadata)){
           print("IndicationOn feature in Metadata indicates whether sample received stimulation or not ")
           Metadata$IndicationOn = 1
           Metadata$IndicationOn[is.na(Metadata$Indication_induction_time_hrs)] = 0
        }else if("stim" %in% colnames(Metadata)){
            Metadata$IndicationOn = 1
            Metadata$IndicationOn[grep("^No",Metadata$stim)] = 0
        }else{

            print("IndicationOn neither Stimulant_used was found within the metadata")
        }

        if("Treatment_conc_uM" %in% colnames(Metadata)){

            # As levels for PCA
            Metadata$Treatment_Conc_uM = factor(Metadata$Treatment_conc_uM, levels = sort(unique(as.numeric(Metadata$Treatment_conc_uM))))
            # As numeric for correlation plots to see them with colors in order
            Metadata$Treatment_conc_uM = as.numeric(Metadata$Treatment_conc_uM)
        }else if("Treatment_concentration_mg_Kg" %in% colnames(Metadata)){
            Metadata$Treatment_Conc_mg = factor(Metadata$Treatment_concentration_mg_Kg, levels = sort(unique(as.numeric(Metadata$Treatment_concentration_mg_Kg))))
            Metadata$Treatment_conc_mg = as.numeric(Metadata$Treatment_conc_mg)
        }else{

            print("Treatment concentration column not found")

        }


        if("RNA.Sample::RQN"  %in% colnames(Metadata)){

            Metadata$RQN = Metadata$`RNA.Sample::RQN`
        }else if ("RQN"  %in% colnames(Metadata)){

            print("RQN column already in Metadata")
        }else{

            print("RQN column not found")

        }

        Metadata$Rna_collection_time_hrs = factor(Metadata$Rna_collection_time_hrs, levels = sort(unique(as.numeric(Metadata$Rna_collection_time_hrs))))


        # ALL SAMPLES

        PCA = Plot_PC_Invivo_PerMol(t(NormCounts), Metadata, Color_gg = 'Treatment', Title = 'PCA all samples')

        # pdf(file= paste(Output_file_path, Name_folder,Name_folder_PCA, "PCA_All_Treatment.png", sep = '/'))
        # print(PCA[[1]])
        # print(PCA[[2]])
        # dev.off()

        P1 = plot_grid(PCA[[1]], PCA[[2]],  ncol=1, align='v')

        ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, "PCA_All_Treatment.pdf", sep = '/'), plot=P1, width = 8,height = 8, units = 'in')
        ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, "PCA_All_Treatment.jpeg", sep = '/'), plot=P1, width = 8,height = 8, units = 'in')

        PCA2 = Plot_PC_Invivo_PerMol(t(NormCounts), Metadata, Color_gg = 'Stimulant_used', Title = 'PCA all samples')

        P2 = plot_grid(PCA2[[1]], PCA2[[2]],  ncol=1, align='v')

        ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, "PCA_All_Stim.pdf", sep = '/'), plot=P2, width = 8,height = 8, units = 'in')
        ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, "PCA_All_Stim.jpeg", sep = '/'), plot=P2, width = 8,height = 8, units = 'in')


        PCA3 = Plot_PC_Invivo_PerMol(t(NormCounts), Metadata, Color_gg = 'RQN', Title = 'PCA all samples')

        P3 = plot_grid(PCA3[[1]], PCA3[[2]],  ncol=1, align='v')

        ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, "PCA_All_RQN.pdf", sep = '/'), plot=P3, width = 8,height = 8, units = 'in')
        ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, "PCA_All_RQN.jpeg", sep = '/'), plot=P3, width = 8,height = 8, units = 'in')

        PCA4 = Plot_PC_Invivo_PerMol(t(NormCounts), Metadata, Color_gg = 'Outliers', Title = 'PCA all samples')

        P4 = plot_grid(PCA4[[1]], PCA4[[2]],  ncol=1, align='v')

        ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, "PCA_All_Outliers.pdf", sep = '/'), plot=P4, width = 8,height = 8, units = 'in')
        ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, "PCA_All_Outliers.jpeg", sep = '/'), plot=P4, width = 8,height = 8, units = 'in')

        PCA5 = Plot_PC_Invivo_PerMol(t(NormCounts), Metadata, Color_gg = 'Plate.id', Title = 'PCA all samples')

        P5 = plot_grid(PCA5[[1]], PCA5[[2]],  ncol=1, align='v')

        ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, "PCA_All_Plate_id.pdf", sep = '/'), plot=P5, width = 8,height = 8, units = 'in')
        ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, "PCA_All_Plate_id.jpeg", sep = '/'), plot=P5, width = 8,height = 8, units = 'in')



        # By time point
        Times = unique(Metadata$timeColl)

        for(Time in Times){

            Metadata_Time = Metadata %>% dplyr::filter(timeColl == Time)
            NormCounts_Time = NormCounts[ ,match(Metadata_Time$Sample_name,colnames(NormCounts) ) ]

            # Get rid of columns with variance equal to zero
            vars = apply(NormCounts_Time, 1, var)
            NormCounts_Time =NormCounts_Time[vars > 0, ]

            PCA1 = Plot_PC_Invivo_PerMol(t(NormCounts_Time), Metadata_Time, Color_gg = 'Treatment', Title = paste('PCA: Time point',Time, sep = ' ' ))

            PT1 = plot_grid(PCA1[[1]], PCA1[[2]],  ncol=1, align='v')

            ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATime", Time, "Treatment",".pdf", sep = '_'), sep = '/'), plot=PT1, width = 8,height = 8, units = 'in')
            ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATime", Time, "Treatment",".jpeg", sep = '_'), sep = '/'), plot=PT1, width = 8,height = 8, units = 'in')


            PCA2 = Plot_PC_Invivo_PerMol(t(NormCounts_Time), Metadata_Time, Color_gg = 'Stimulant_used', Title = paste('PCA: Time point',Time, sep = ' ' ))

            PT2 = plot_grid(PCA2[[1]], PCA2[[2]],  ncol=1, align='v')

            ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATime", Time, "Stim",".pdf", sep = '_'), sep = '/'), plot=PT2, width = 8,height = 8, units = 'in')
            ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATime", Time, "Stim",".jpeg", sep = '_'), sep = '/'), plot=PT2, width = 8,height = 8, units = 'in')


            PCA3 = Plot_PC_Invivo_PerMol(t(NormCounts_Time), Metadata_Time, Color_gg = 'RQN', Title = paste('PCA: Time point',Time, sep = ' ' ))

            PT3 = plot_grid(PCA3[[1]], PCA3[[2]],  ncol=1, align='v')

            ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATime", Time, "RQN",".pdf", sep = '_'), sep = '/'), plot=PT3, width = 8,height = 8, units = 'in')
            ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATime", Time, "RQN",".jpeg", sep = '_'), sep = '/'), plot=PT3, width = 8,height = 8, units = 'in')

            PCA4 = Plot_PC_Invivo_PerMol(t(NormCounts_Time), Metadata_Time, Color_gg = 'Outliers', Title = paste('PCA: Time point',Time, sep = ' ' ))

            PT4 = plot_grid(PCA4[[1]], PCA4[[2]],  ncol=1, align='v')

            ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATime", Time, 'Outliers',".pdf", sep = '_'), sep = '/'), plot=PT4, width = 8,height = 8, units = 'in')
            ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATime", Time, 'Outliers',".jpeg", sep = '_'), sep = '/'), plot=PT4, width = 8,height = 8, units = 'in')

            PCA5 = Plot_PC_Invivo_PerMol(t(NormCounts_Time), Metadata_Time, Color_gg = 'Plate.id', Title = 'PCA all samples')

            PT5 = plot_grid(PCA5[[1]], PCA5[[2]],  ncol=1, align='v')

            ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATime", Time, 'Plate_id',".pdf", sep = '_'), sep = '/'), plot=PT5, width = 8,height = 8, units = 'in')
            ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATime", Time, 'Plate_id',".jpeg", sep = '_'), sep = '/'), plot=PT5, width = 8,height = 8, units = 'in')


            # By molecule
            Mols = unique(Metadata_Time$Treatment)

            for(Mol in Mols){

                Metadata_Time_mol = Metadata_Time %>% dplyr::filter(Treatment == Mol)
                NormCounts_Time_mol = NormCounts_Time[ ,match(Metadata_Time_mol$Sample_name,colnames(NormCounts_Time) ) ]

                # Get rid of columns with variance equal to zero
                vars = apply(NormCounts_Time_mol, 1, var)
                NormCounts_Time_mol =NormCounts_Time_mol[vars > 0, ]

                PCA_M1 = Plot_PC_Invivo_PerMol(t(NormCounts_Time_mol), Metadata_Time_mol, Color_gg = 'Stimulant_used', Title = paste('PCA: Time point',Time,', Mol',Mol, sep = ' ' ))

                PTM1 = plot_grid(PCA_M1[[1]], PCA_M1[[2]],  ncol=1, align='v')

                ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATimeMol", Time,Mol, "Stim",".pdf", sep = '_'), sep = '/'), plot=PTM1, width = 8,height = 8, units = 'in')
                ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATimeMol", Time,Mol, "Stim",".jpeg", sep = '_'), sep = '/'), plot=PTM1, width = 8,height = 8, units = 'in')

                PCA_M2 = Plot_PC_Invivo_PerMol(t(NormCounts_Time_mol), Metadata_Time_mol, Color_gg = 'Treatment_Conc_uM', Title = paste('PCA: Time point',Time,', Mol',Mol, sep = ' ' ))

                PTM2 = plot_grid(PCA_M2[[1]], PCA_M2[[2]],  ncol=1, align='v')

                ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATimeMol", Time,Mol, "Conc",".pdf", sep = '_'), sep = '/'), plot=PTM2, width = 8,height = 8, units = 'in')
                ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATimeMol", Time,Mol, "Conc",".jpeg", sep = '_'), sep = '/'), plot=PTM2, width = 8,height = 8, units = 'in')

                PCA_M3 = Plot_PC_Invivo_PerMol(t(NormCounts_Time_mol), Metadata_Time_mol, Color_gg = 'RQN', Title = paste('PCA: Time point',Time,', Mol',Mol, sep = ' ' ))

                PTM3 = plot_grid(PCA_M3[[1]], PCA_M3[[2]],  ncol=1, align='v')

                ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATimeMol", Time,Mol, 'RQN',".pdf", sep = '_'), sep = '/'), plot=PTM3, width = 8,height = 8, units = 'in')
                ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATimeMol", Time,Mol, 'RQN',".jpeg", sep = '_'), sep = '/'), plot=PTM3, width = 8,height = 8, units = 'in')


                PCA_M4 = Plot_PC_Invivo_PerMol(t(NormCounts_Time_mol), Metadata_Time_mol, Color_gg = 'Outliers', Title = paste('PCA: Time point',Time,', Mol',Mol, sep = ' ' ))

                PTM4 = plot_grid(PCA_M4[[1]], PCA_M4[[2]],  ncol=1, align='v')

                ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATimeMol", Time,Mol, 'Outliers',".pdf", sep = '_'), sep = '/'), plot=PTM4, width = 8,height = 8, units = 'in')
                ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATimeMol", Time,Mol, 'Outliers',".jpeg", sep = '_'), sep = '/'), plot=PTM4, width = 8,height = 8, units = 'in')

                PCA_M5 = Plot_PC_Invivo_PerMol(t(NormCounts_Time_mol), Metadata_Time_mol, Color_gg = 'Plate.id', Title = paste('PCA: Time point',Time,', Mol',Mol, sep = ' ' ))

                PTM5 = plot_grid(PCA_M5[[1]], PCA_M5[[2]],  ncol=1, align='v')

                ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATimeMol", Time,Mol, 'Plate_id',".pdf", sep = '_'), sep = '/'), plot=PTM5, width = 8,height = 8, units = 'in')
                ggsave(filename=paste(Output_file_path, Name_folder,Name_folder_PCA, paste("PCATimeMol", Time,Mol, 'Plate_id',".jpeg", sep = '_'), sep = '/'), plot=PTM5, width = 8,height = 8, units = 'in')


            }



        }


        # CORRELATION PLOTS

        annotation_col = as.data.frame(Metadata %>% dplyr::select(Rna_collection_time_hrs, Stimulant_used, Treatment, Treatment_conc_uM, Sample_name))
        rownames(annotation_col) = annotation_col$Sample_name
        annotation_col$Sample_name = NULL

        annotation_row = as.data.frame(Metadata %>% dplyr::select(Sample_name,  RQN, Outliers))
        rownames(annotation_row) = annotation_row$Sample_name
        annotation_row$Sample_name = NULL

        corr_sp <- cor(NormCounts, method = "spearman")
        paletteLength <- 50
        myColor <- viridis::viridis(paletteLength)
        ann_colors = list(
            Outliers = c('No' = "navy", 'Yes' = "firebrick")
        )

        jpeg(file = paste(Output_file_path, Name_folder,Name_folder_PCA, "CorrSp_All.jpeg", sep = '/'), width = 10, height = 8, units = "in", res = 300)

            P = pheatmap(corr_sp, annotation_col = annotation_col,annotation_row = annotation_row, show_colnames = F, show_rownames = F,
                     clustering_method = "ward.D2",fontsize_row =5, color = rev(myColor), annotation_colors = ann_colors)
            print(P)

        dev.off()

        corr_pr <- cor(NormCounts, method = "pearson")

        jpeg(file = paste(Output_file_path, Name_folder,Name_folder_PCA, "CorrPr_All.jpeg", sep = '/'), width = 10, height = 8, units = "in", res = 300)

           P2 = pheatmap(corr_pr, annotation_col = annotation_col,annotation_row = annotation_row, show_colnames = F, show_rownames = F,
                     clustering_method = "ward.D2",fontsize_row =5, color = rev(myColor), annotation_colors = ann_colors)
           print(P2)

        dev.off()


        # Show relationship between a handful samples

        outliers_left = outliers_cust[outliers_cust %in% colnames(NormCounts)]
        if(length(outliers_left) != 0 ){
            sam = sample(colnames(NormCounts)[!colnames(NormCounts) %in% outliers_cust], 3)
            sam = c(sam, outliers_left[1:min(3, length(outliers_left))])

        }else{

            sam = sample(colnames(NormCounts)[!colnames(NormCounts) %in% outliers_cust], 6)
        }


        jpeg(file = paste(Output_file_path, Name_folder,Name_folder_PCA, "PairSp_RandSam.jpeg", sep = '/'), width = 10, height = 8, units = "in", res = 300)

           P3 = ggpairs(NormCounts[,sam], upper = list(continuous = wrap(ggally_cor, method = "spearman"))) + theme(strip.text.x = element_text(size = 5),
                                                                                                               strip.text.y = element_text(size = 5))
           print(P3)

        dev.off()

        # By time point

        Times = unique(Metadata$timeColl)
        for(Time in Times){

            Metadata_Time = Metadata %>% dplyr::filter(timeColl == Time)
            NormCounts_Time = NormCounts[ ,match(Metadata_Time$Sample_name,colnames(NormCounts) ) ]

            annotation_col = as.data.frame(Metadata_Time %>% dplyr::select( Stimulant_used, Treatment, Treatment_conc_uM, Sample_name))
            rownames(annotation_col) = annotation_col$Sample_name
            annotation_col$Sample_name = NULL

            annotation_row = as.data.frame(Metadata_Time %>% dplyr::select(Sample_name,  RQN, Outliers))
            rownames(annotation_row) = annotation_row$Sample_name
            annotation_row$Sample_name = NULL

            corr_sp_time <- cor(NormCounts_Time, method = "spearman")
            paletteLength <- 50
            myColor <- viridis::viridis(paletteLength)
            ann_colors = list(
                Outliers = c('No' = "navy", 'Yes' = "firebrick")
            )

            jpeg(file = paste(Output_file_path, Name_folder,Name_folder_PCA, paste("CorrSpTime",Time,".jpeg",sep = '_'), sep = '/'), width = 10, height = 8, units = "in", res = 300)

                  P = pheatmap(corr_sp_time, annotation_col = annotation_col,annotation_row = annotation_row, show_colnames = F, show_rownames = F,
                         clustering_method = "ward.D2",fontsize_row =5, color = rev(myColor), annotation_colors = ann_colors)
                  print(P)

            dev.off()

            corr_pr_time <- cor(NormCounts_Time, method = "pearson")

            jpeg(file = paste(Output_file_path, Name_folder,Name_folder_PCA, paste("CorrPrTime",Time,".jpeg",sep = '_'), sep = '/'), width = 10, height = 8, units = "in", res = 300)

                 P2 = pheatmap(corr_pr_time, annotation_col = annotation_col,annotation_row = annotation_row, show_colnames = F, show_rownames = F,
                          clustering_method = "ward.D2",fontsize_row =5, color = rev(myColor), annotation_colors = ann_colors)
            print(P2)

            dev.off()



        }




    }

    setwd(Output_file_path)

    return()
}
