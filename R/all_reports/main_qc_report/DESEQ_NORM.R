# DESeq normalization and DEG
DESEQ_NORM = function(Output_file_path, QCNORM = "PRE_QCNORM", outliers_path, CoarseConditions, METHOD_NORM = 'Standard', VST_FILTER = "VST_ON",  SVD_FILTER = 'SVD_OFF', MEM_MB = 1600, Batch_variable = 'Plate.id', DataForPipeline ){

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
    library(limma)
    library(edgeR)



    #source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/Functions_Invivo.R")

    clusterTemp = file.path(DataForPipeline, "slurmMqBoris.tmpl")

    # Import data
    setwd(Output_file_path)
    Data_file = "./PRE_FILTERING/Prefilter_Data.xlsx"
    Count = read.xlsx(xlsxFile = Data_file, sheet = "Count", rowNames= TRUE)
    Metadata = read.xlsx(xlsxFile = Data_file, sheet = "Metadata")

    #Metadata$Treatment_conc = as.numeric(Metadata$Treatment_conc)


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
    # if(GENE_ENSEM == TRUE){
    #     if(ORGANISM == 'Human'){
    #     # Convert from Ensemble to symbol
    #     geneSymbol <- mapIds(org.Hs.eg.db,keys=rownames(Count),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
    #     Count$genesymbol = geneSymbol[match(rownames(Count), names(geneSymbol))]
    #     Count = Count[!is.na(Count$genesymbol),]
    #
    #     Count <- aggregate(. ~ genesymbol, data = Count, FUN = sum)
    #     rownames(Count) = Count$genesymbol
    #     Count$genesymbol = NULL
    #   }
    # }


    # Normalization
    if (METHOD_NORM == 'Standard'){
        print('Running DESEQ Standard normalization...')
        ##create a DESeq object

        Design = "~ CoarseCondition"
        # To be modified if more than one factor is desired to be adjusted for
        if(Batch_variable != "" ){

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

        #Save RData object
        save(deseqObj, NormCounts,Metadata, file = 'DEG_Norm.RData')

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
            if(Batch_variable != "" ){

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

            #Save RData object
            save(deseqObj, NormCounts,Metadata, file = 'DEG_Norm.RData')

    }else if(METHOD_NORM == "limma_voom"){

        # To be modified if more than one factor is desired to be adjusted for
        design <- model.matrix(~ CoarseCondition, data = Metadata)
        dge <- DGEList(counts = Count)
        dge <- calcNormFactors(dge, method="TMM")
        v <- voom(dge, design)
        NormCounts = v$E

        save(dge, NormCounts,Metadata, file = 'DEG_Norm.RData')

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



    # Save excel file
    wb <- createWorkbook()
    addWorksheet(wb, "NormCounts")
    addWorksheet(wb, "Metadata")
    writeData(wb, sheet = "NormCounts", x = NormCounts, rowNames = TRUE)
    writeData(wb, sheet = "Metadata", x = Metadata)
    saveWorkbook(wb, "Norm_Data.xlsx", overwrite = TRUE)


    setwd(Output_file_path)

    return()
}
