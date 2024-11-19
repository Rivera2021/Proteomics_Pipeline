
# Modify metadata to have a colum name that matches sample names in count matrix
Metadata_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/In_vitro/Watchmaker 25_50_Metadata_Standarized.xlsx"
Metadata = read_xlsx(Metadata_path)

Input_matrix_path = "/fsx/home/crivera/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/In_vitro/salmon.merged.gene_counts.tsv"
Count = read.table(Input_matrix_path,  sep = '\t', header = TRUE, check.names=FALSE)
Original_colnames = colnames(Count)[3:length(colnames(Count))]
Wellcolnames =   unlist(lapply(colnames(Count), function(x){str_split(x,'_')[[1]][2]}))
Wellcolnames = Wellcolnames[3: length(Wellcolnames)]
Df = data.frame(Well = Wellcolnames)
Df$Original = Original_colnames

Metadata$Scinamic = Metadata$`ELN_ ID`
Metadata_temp = Metadata %>% dplyr::filter(Scinamic == "A-2024-0142")
Metadata_temp$OriginalName = Df$Original[match(Metadata_temp$Library_Plate_Well, Df$Well)]

write.xlsx(Metadata_temp, '~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/In_vitro/Metadata_mod.xlsx')
