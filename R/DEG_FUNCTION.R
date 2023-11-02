# DEG per time point and Pathway Enrichment
DEG_FUNCTION = function(Output_file_path, List_contrasts_Path, DEG_Method = 'DESeq',MH_Method = 'BH', AlphaHC = 0.1,  padjval = 0.2 ,LogFoldThrs_VolPlot = 1){
    library(ggplot2)
    library(ggrepel)
    library(DESeq2)
    library(ddpca)
    library(foreach)
    library(BioMark)

    source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/Functions_Invivo.R")
    # Intro message
    print("Getting DEG...")

    # Import data
    setwd(Output_file_path)
    Data_file = "./DESEQ_NORM/DESeq_Norm.RData"
    load(Data_file)

    List_contrasts = read.xlsx(xlsxFile = List_contrasts_Path)

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
        wb <- createWorkbook()

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

        for(i in 1:nrow(List_contrasts_mol)){
            print(List_contrasts_mol$treat[i])
            Name_contr = paste(List_contrasts_mol$treat[i],'VS',List_contrasts_mol$untreat[i],sep = '_')
            print(Name_contr)
            # Filter necessary data for comparison from Metadata
            Metadata_temp = Metadata %>% filter(CoarseCondition %in% c(List_contrasts_mol$treat[i],List_contrasts_mol$untreat[i]))
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
                ModelInfo_sig = ModelInfo %>% filter(pvajSig == 1)

            }else if(MH_Method == 'BH'){
                if(DEG_Method != 'DESeq'){
                    ModelInfo$padj =p.adjust(as.numeric(ModelInfo$p.value), method='BH')
                    ModelInfo_sig = ModelInfo %>% filter(padj < padjval)
                }else if(DEG_Method == 'DESeq'){

                    ModelInfo_sig = ModelInfo %>% filter(padj < padjval)
                    print(nrow(ModelInfo_sig))
                }

                # Plot volcano plot
                VP = Volcano_plot_padj(res_Vol = ModelInfo, padj_thr = padjval, logFC_thr = LogFoldThrs_VolPlot, Title = Name_contr , maxOver = 20,xLabel = 'Log2Fold',yLabel = '-Log10(P-adj)')

                png(filename=paste("VolcanoPlot", Name_contr,".png", sep = '_'), width = 1000, height = 1000, res=100)
                print(VP)
                dev.off()

            }


            DiffGenes[[Name_contr]] = ModelInfo


            addWorksheet(wb, Name_contr)
            writeData(wb, sheet = Name_contr, x = ModelInfo_sig)


        }

        save(DiffGenes, file = paste('DEG', drug,'.RData' , sep = '_'))
        saveWorkbook(wb, paste('DEG', drug,'.xlsx', sep = '_'), overwrite = TRUE)

        setwd(paste(Output_file_path,Name_folder, sep = '/'))


    }


}
