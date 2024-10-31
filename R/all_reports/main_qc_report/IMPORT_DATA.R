

# Reformat raw counts and metadata
IMPORT_DATA = function(Output_file_path, Metadata_path, CountM_path, Metrics_path, TPM_path, Organism_type, SampleName_path){

    library("openxlsx")
    library("readxl")
    library(tidyverse)

    dir.create(Output_file_path)
    setwd(Output_file_path)

    # Save metadata and count formatted data frames
    dir.create("IMPORT_DATA")
    setwd("IMPORT_DATA")

    if(file.exists(Metadata_path)){
    Metadata = read_xlsx(Metadata_path)
    Metadata$scinamicNum = Metadata$`Biosample::Experiment Number`


    # Modify Treatment notation for convenience
    Metadata$Treatment = Metadata$`Biosample::Compound Batch::ID`
    Metadata$Stimulant_treatment_time_hrs= Metadata$`Biosample::Stimulant Treatment Time (h)`
    Metadata$Replicate = Metadata$`Biosample::Replicate`

    # Replace NA by zero in Treatment concentration. Order concentrations
    Metadata$Treatment_conc = Metadata$`Biosample::Compound Concentration`
    if(Organism_type == "In_vivo"){
        Metadata$Stimulant_used = ifelse(Metadata$`Biosample::Compound Batch::ID`!="NONE", "Stim", "No_Stim")
        Metadata$Treatment_time_hrs = Metadata$`Animal Tissue::Time of Collection (h)`
        Metadata$RQN = Metadata$`RNA Sample::RQN`
        Metadata$RNA_Conc = Metadata$`RNA Sample::Concentration (ng/uL)`

    }else if(Organism_type == "In_vitro"){

        Metadata$Stimulant_used = Metadata$`Biosample::Stimulant`
        Metadata$Treatment_time_hrs = Metadata$`Biosample::Compound Treatment Time (h)`
        Metadata$RQN = Metadata$RQN
        Metadata$RNA_Conc = Metadata$`Concentration (ng/uL)`
        Metadata$Cell_line = Metadata$`Biosample::Cellline::ID`
    }else{
        stop("Please provide a valid Organism")
    }


    SampleNames <- read.csv(SampleName_path, header = FALSE)
    SampleNames = SampleNames$V1
    if(any(!SampleNames %in%colnames(Metadata))){
        stop("Choose valid variables for sample names")
    }else{

        Metadata$Sample_name = Metadata[[SampleNames[1]]]

        for(j in 2:length(SampleNames)){
            Metadata$Sample_name = paste(Metadata[["Sample_name"]],Metadata[[SampleNames[j]]], sep = "_")
        }
    }
    }else{
        stop("Metadata file does not exist")
    }

    # Import Count Matrix
    if(file.exists(CountM_path)){
    Count = read.table(CountM_path,  sep = '\t', header = TRUE, check.names=FALSE)
    Count$gene_id = NULL
    # Aggregate reads for duplicated genes
    Count <- aggregate(. ~ gene_name, data = Count, FUN = sum)
    Count = Count %>% column_to_rownames("gene_name")

    if(any(!colnames(Count) %in% Metadata$ID)){

        print("There are samples in the count matrix that do not match any sample in the metadata")
    }
    if(any(!Metadata$ID %in% colnames(Count))){

        print("There are samples in the metadata that do not match any column in the count matrix")
    }

    Metadata = Metadata %>% dplyr::filter(ID  %in% colnames(Count))
    Count = Count[,match(Metadata$ID, colnames(Count))]
    }else{
        stop("Path for raw count matrix does not exist")
    }


    # Read TPM matrix
    if(file.exists(TPM_path)){
    TPM = read.table(TPM_path,  sep = '\t', header = TRUE, check.names=FALSE)
    TPM$gene_id = NULL
    TPM <- aggregate(. ~ gene_name, data = TPM1, FUN = sum)
    TPM = TPM %>% column_to_rownames("gene_name")
    TPM = TPM[,match(Metadata$ID, colnames(TPM))]

    write.table(TPM, file='./TPM.tsv', row.names = FALSE,quote=FALSE, sep='\t')

    }else{
        print("TPM file does not exist")
    }
    # Read MultiQC files
    if (file.exists(Metrics_path)){
        multiqc_metrics <- rio::import(Metrics_path)%>% janitor::clean_names() %>% as_tibble()
        multiqc_metrics = multiqc_metrics[!is.na(multiqc_metrics$custom_content_biotype_counts_mqc_generalstats_custom_content_biotype_counts_percent_r_rna),]

        write.table(multiqc_metrics, file='./MultiQC.tsv', row.names = FALSE,quote=FALSE, sep='\t')

    }else{
        print("Path for file with QC metrics does not exist")
    }

    # Save metadata and count formatted data frames

    wb <- createWorkbook()
    addWorksheet(wb, "Count")
    addWorksheet(wb, "Metadata")

    writeData(wb, sheet = "Count", x = Count)
    writeData(wb, sheet = "Metadata", x = Metadata)

    saveWorkbook(wb, "Import_Data.xlsx", overwrite = TRUE)
    setwd(Output_file_path)




}
