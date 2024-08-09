# DEG per time point and Pathway Enrichment
DEG_FUNCTION = function(Output_file_path, outliers_path, DEG_QCNORM,Control_treatment  ,  List_contrasts_Path, DEG_Method = 'DESeq',MH_Method = 'BH', AlphaHC = 0.1,  padjval = 0.2 ,LogFoldThrs_VolPlot = 1){
  library(ggplot2)
  library(ggrepel)
  library(DESeq2)
  library(ddpca)
  library(foreach)
  library(BioMark)
  library(stringr)
  library("openxlsx")
  library(dplyr)
  library(org.Hs.eg.db)

  source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/Functions_Invivo.R")
  # Intro message
  print("Getting DEG...")

  # Import data

  setwd(Output_file_path)
  Data_file = "./PRE_FILTERING/Prefilter_Data.xlsx"
  Count = read.xlsx(xlsxFile = Data_file, sheet = "Count", rowNames= TRUE)
  Metadata = read.xlsx(xlsxFile = Data_file, sheet = "Metadata")

  # Convert from Ensemble to rownames
  geneSymbol <- mapIds(org.Hs.eg.db,keys=rownames(Count),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
  Count$genesymbol = geneSymbol[match(rownames(Count), names(geneSymbol))]
  Count = Count[!is.na(Count$genesymbol),]

  Count <- aggregate(. ~ genesymbol, data = Count, FUN = sum)
  rownames(Count) = Count$genesymbol
  Count$genesymbol = NULL

  # Read customize outliers
  if(file.exists(outliers_path)){

      outliers_cust <- read.csv(outliers_path, header = FALSE)
      outliers_cust = outliers_cust[[1]]
  }else{

      print("Path to costumized outliers does not exist or is null")
      outliers_cust = c()
  }

  # Evaluate condition depending on whether the function is call before or after QC_POSTNORNMALIZATION
  if(DEG_QCNORM == "DEG_POST_QCNORM"){
      outliers = outliers_cust

      # Create Folder

      Name_folder =  "DEG_POST_QCNORM"
      dir.create(Name_folder)
      setwd(Name_folder)

  }else if(DEG_QCNORM == "DEG_PRE_QCNORM"){

      Name_folder = "DEG_PRE_QCNORM"
      dir.create(Name_folder)
      setwd(Name_folder)

      outliers = outliers_cust

  }else{

      stop("Select a valid DEG_QCNORM option")
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


  print('Running DESEQ Standard normalization...')
  # Use this more elegant way in the future
  #grouped_list <- split(Metadata, interaction(Metadata$Treatment, Metadata$Cell_line, Metadata$Stimulant_used, Metadata$timeColl, Metadata$Treatment_conc_uM))

  drugs = unique(Metadata$Treatment)
  drugs = drugs[drugs != Control_treatment]
  for(mol in drugs){

      Metadata_mol = Metadata %>% filter(Treatment == mol)
      Count_mol = Count[ ,match(Metadata_mol$Sample_name,colnames(Count))]

      cell_lines = unique(Metadata_mol$Cell_line)
      for(cellline in cell_lines){

          Metadata_mol_cell = Metadata_mol %>% filter(Cell_line == cellline)
          Count_mol_cell = Count_mol[ ,match(Metadata_mol_cell$Sample_name,colnames(Count_mol))]


          stims = unique(Metadata_mol_cell$Stimulant_used)

          for(stim in stims){

             Metadata_mol_cell_stim = Metadata_mol_cell %>% filter(Stimulant_used == stim)
             Count_mol_cell_stim = Count_mol_cell[ ,match(Metadata_mol_cell_stim$Sample_name,colnames(Count_mol_cell))]

             Times = unique(Metadata_mol_cell_stim$timeColl)

             for(time in Times){

                  Metadata_mol_cell_stim_time = Metadata_mol_cell_stim %>% filter(timeColl == time)
                  Count_mol_cell_stim_time  = Count_mol_cell_stim[ ,match(Metadata_mol_cell_stim_time$Sample_name,colnames(Count_mol_cell_stim))]

                  doses = unique(Metadata_mol_cell_stim_time$Treatment_conc_uM)

                  for(dose in doses){

                   Metadata_mol_cell_stim_time_dose = Metadata_mol_cell_stim_time %>% filter(Treatment_conc_uM == dose)
                   Count_mol_cell_stim_time_dose = Count_mol_cell_stim_time[ ,match(Metadata_mol_cell_stim_time_dose$Sample_name,colnames(Count_mol_cell_stim_time))]


                   Metadata_comp = Metadata %>% filter(Treatment %in% c(mol,Control_treatment) & Cell_line == cellline & Stimulant_used == stim & timeColl == time & Treatment_conc_uM %in% c(dose,0 ))
                   Count_comp  = Count[ ,match(Metadata_comp$Sample_name,colnames(Count))]

                   # DESeq DEG
                   dds <- DESeqDataSetFromMatrix(countData = Count_comp, colData = Metadata_comp, design = ~Treatment)

                   deseqObj = DESeq(dds,
                                    parallel = TRUE,
                                    fitType = "parametric",
                                    BPPARAM=MulticoreParam(detectCores() - 2))

                   res = results(deseqObj,contrast = c("Treatment",mol, Control_treatment ) )
                   ModelInfo = data.frame(res)
                   ModelInfo$genes = rownames(ModelInfo)

                   file_name = paste('mol', mol, 'cellln', cellline, 'stim', stim, 'time', time,'dose', dose, '.xlsx', sep ='_' )
                   write.xlsx(ModelInfo, file = file_name)
                   print(file_name)
                   ModelInfo_sig = ModelInfo %>% filter(padj < 0.2)

                   print(nrow(ModelInfo_sig))


                 }

             }

          }
      }


  }

  # DEG for stimulation

  Name_folder = "STIMULATION"
  dir.create(Name_folder)
  setwd(Name_folder)


  NOStim = 'NONE'
  stims = c("R848", "TNFa")
  #stims = stims[stims!= NOStim]
  Treatments = c('DMSO', 'NONE')

  Diff_Stim = list()
  wb <- createWorkbook()

  for(stim in stims){


      print(stim)
      Metadata_stim = Metadata %>% filter(Stimulant_used == stim)
      Count_stim = Count[,match(Metadata_stim$Sample_name,colnames(Count) )]

      celllines = unique(Metadata_stim$Cell_line)

      for(cellline in celllines){

          print(cellline)
          Metadata_stim_cell = Metadata_stim %>% filter(Cell_line == cellline)
          Count_stim_cell = Count_stim[,match(Metadata_stim_cell$Sample_name,colnames(Count_stim) )]

          Times = unique( Metadata_stim_cell$timeColl)
          for(Time in Times){

              Metadata_stim_cell_Time = Metadata_stim_cell %>% filter(timeColl == Time)
              Count_stim_cell_Time = Count_stim_cell[,match(Metadata_stim_cell_Time$Sample_name,colnames(Count_stim_cell))]


              for(mol in Treatments){
                  Metadata_comp = Metadata %>% filter(Stimulant_used %in% c(stim, NOStim) & Cell_line == cellline & timeColl == Time & Treatment == mol)
                  Count_comp  = Count[ ,match(Metadata_comp$Sample_name,colnames(Count))]


                  Count_comp_filt = Count_comp[rowSums(Count_comp) > ceiling(ncol(Count_comp)*0.2), ]
                  # DESeq DEG

                  dds <- DESeqDataSetFromMatrix(countData = Count_comp_filt, colData = Metadata_comp, design = ~Stimulant_used)

                  deseqObj = DESeq(dds,
                                   parallel = TRUE,
                                   fitType = "parametric",
                                   BPPARAM=MulticoreParam(detectCores() - 2))

                  resT = results(deseqObj,contrast = c("Stimulant_used", stim, NOStim ) )
                  ModelInfoT = data.frame(resT)
                  ModelInfoT$Genes = rownames(ModelInfoT)

                  Name = paste( stim, cellline, Time,mol, sep ='_' )
                  Diff_Stim[[Name]] = ModelInfoT %>% filter(padj < 0.2)

                  ModelInfoSig = ModelInfoT %>% filter(padj < 0.2)
                  addWorksheet(wb, Name)
                  writeData(wb, sheet = Name, x = ModelInfoSig)



              }

            }

          }

      }


   saveWorkbook(wb, paste('DEG_DESEQM1_stim_All.xlsx', sep = '_'), overwrite = TRUE)
   save(Diff_Stim, file = paste('DEG_DESEQM1_stim_All.RData' , sep = '_'))

  }







  ##create a DESeq object
  dds <- DESeqDataSetFromMatrix(countData = Count, colData = Metadata, design = ~ CoarseCondition)
  startTime <- Sys.time()
  deseqObj = DESeq(dds,
                   parallel = TRUE,
                   fitType = "parametric",
                   BPPARAM=MulticoreParam(detectCores() - 2))














Data_file = "./DESEQ_NORM_QCNORM/DESeq_Norm.RData"
  load(Data_file)
  #
  #print("Only contrasts allowed after removing outliers will be calculated")
  #List_contrasts = read.xlsx(xlsxFile = List_contrasts_Path)
  #List_contrasts = List_contrasts %>% filter(treat %in% Metadata$CoarseCondition & untreat %in% Metadata$CoarseCondition)

  #drugsInContrasts = unique(sapply(List_contrasts$treat, function(x){str_split(x,'_')[[1]][1]}))
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

      }else{

        stop("Select a valid Multiple hypothesis method")
      }


      DiffGenes[[Name_contr]] = ModelInfo


      addWorksheet(wb, Name_contr)
      writeData(wb, sheet = Name_contr, x = ModelInfo_sig)


    }

    save(DiffGenes, file = paste('DEG', drug,'.RData' , sep = '_'))
    saveWorkbook(wb, paste('DEG', drug,'.xlsx', sep = '_'), overwrite = TRUE)

    setwd(paste(Output_file_path,Name_folder, sep = '/'))


  }
  setwd(Output_file_path)
  return()
}
