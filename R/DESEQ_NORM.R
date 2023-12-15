# DESeq normalization and DEG
DESEQ_NORM = function(Output_file_path, QCNORM = "PRE_QCNORM", outliers_path, CoarseConditions, METHOD_NORM = 'Standard', VST_FILTER = "VST_ON",  SVD_FILTER = 'SVD_OFF', PlotPCA = 'PCA_PLOT', Control_Neg_PCA = "", Control_Pos_PCA = ""){

    library(DESeq2)
    #library(readxl)
    library(dplyr)
    library(BiocParallel)
    library(parallel)
    library("PCAtools")
    library(foreach)
    library(doParallel)
    library(openxlsx)
    nCores = detectCores()
    cl <- makeCluster(nCores)
    registerDoParallel(cl)

    source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/Functions_Invivo.R")

    # Import data

    setwd(Output_file_path)
    Data_file = "./PRE_FILTERING/Prefilter_Data.xlsx"
    Count = read.xlsx(xlsxFile = Data_file, sheet = "Count", rowNames= TRUE)
    Metadata = read.xlsx(xlsxFile = Data_file, sheet = "Metadata")

    # Read customize outliers

    if(file.exists(outliers_path)){

        outliers_cust <- read.csv(outliers_path, header = FALSE)
        outliers_cust = outliers_cust[[1]]
    }else{

        print("Path to costumized outliers does not exists or is null")
        outliers_cust = c()
    }

    # Evaluate condition depending on whether the function is call before or after QC_POSTNORNMALIZATION
    if(QCNORM == "POST_QCNORM"){

        # Join outliers from Outlier_file_path and from QC_POSTNORMALIZATION

        if(class(try(read.csv("./QC_POSTNORMALIZATION/Outliers_Selected.csv", header = TRUE))) == "try-error"){
            print("outlier file from QC_NORMALIZATION is empty. Only customized outliers are considered")
            outliers = outliers_cust

        }else{

            outliers_qc = read.csv("./QC_POSTNORMALIZATION/Outliers_Selected.csv", header = TRUE)
            outliers_qc = outliers_qc[[1]]
            outliers = unique(c(outliers_cust, outliers_qc))

        }

        # Create Folder

        Name_folder =  "DESEQ_NORM_QCNORM"
        dir.create(Name_folder)
        setwd(Name_folder)

    }else if(QCNORM == "PRE_QCNORM"){

        Name_folder =  "DESEQ_NORM"
        dir.create(Name_folder)
        setwd(Name_folder)

        outliers = outliers_cust

    }else{

        stop("Select a valid QC_NORM option")
    }

    # Filter for outliers

    if(length(outliers)>0){
        print('Getting rid of outliers')

        Metadata = Metadata %>% filter(!Sample_name %in% outliers)
        Count = Count[ ,match(Metadata$Sample_name,colnames(Count) ) ]
        print(dim(Count))
    }else{
        print('No outliers list provided')
    }

    # DESeq normalization
    #Create coarse Condition feature and useful Metadata
    Metadata$CoarseCondition = Metadata[[CoarseConditions[1]]]

    for(j in 2:length(CoarseConditions)){
        Metadata$CoarseCondition = paste(Metadata[["CoarseCondition"]],Metadata[[CoarseConditions[j]]], sep = "_")
    }


    # Normalization
    if (METHOD_NORM == 'Standard'){
        print('Running Standard normalization...')
        ##create a DESeq object
        dds <- DESeqDataSetFromMatrix(countData = Count, colData = Metadata, design = ~ CoarseCondition)
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
           Metadata$IndicationOn[is.na(Metadata$`Indication._induction_time(hrs)`)] = 0
        }else if("Stimulant_used" %in% colnames(Metadata)){
            Metadata$IndicationOn = 1
            Metadata$IndicationOn[is.na(Metadata$Stimulant_used)] = 0
        }else{

            print("IndicationOn neither Stimulant_used was found within the metadata")
        }

        if("Treatment_conc_uM" %in% colnames(Metadata)){

            Metadata$Treatment_Conc = Metadata$Treatment_conc_uM
        }else if("Treatment_concentration_mg_Kg" %in% colnames(Metadata)){
            Metadata$Treatment_Conc = Metadata$Treatment_concentration_mg_Kg
        }else{

            print("Treatment concentration column not found")
        }

        PCA = Plot_PC_Invivo_V1(t(NormCounts), Metadata, Color_gg = "Treatment", Title = 'PCA all samples')

        pdf(file= "PCA_allSamples.pdf")
        print(PCA[[1]])
        print(PCA[[2]])
        print(PCA[[3]])
        dev.off()

        # Per Mol
        Drugs = unique(Metadata$Treatment)
        Drugs = Drugs[grep('M',Drugs)]

        # Rename time to be able to plot in ggplot

        #Metadata$timeColl = factor(Metadata$`Rna_collection_time(hrs)`, levels = sort(as.numeric(unique(Metadata$`Rna_collection_time(hrs)`))))



        for(drug in Drugs){

            Metadata_mol = Metadata %>% filter(Treatment %in% c(drug, Control_Neg_PCA, Control_Pos_PCA))
            NormCounts_mol = NormCounts[ ,match(Metadata_mol$Sample_name,colnames(NormCounts) ) ]

            # Get rid of columns with variance equal to zero
            vars = apply(NormCounts_mol, 1, var)
            NormCounts_mol =NormCounts_mol[vars > 0, ]


            PCA = Plot_PC_Invivo_PerMol(t(NormCounts_mol), Metadata_mol, Color_gg = 'Treatment', Title = paste('PCA', drug, sep = '_'))

            pdf(file=paste("PCA", drug, ".pdf", sep = '_'))
            print(PCA[[1]])
            print(PCA[[2]])
            print(PCA[[3]])
            print(PCA[[4]])
            print(PCA[[5]])
            print(PCA[[6]])
            dev.off()

            png(filename=paste("PCA", drug,".png", sep = '_'), width = 1000, height = 1000, res=100)
            print(PCA[[1]])
            dev.off()


        }






    }

    setwd(Output_file_path)

    return()
}
