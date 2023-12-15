# Introducing some extra-prefiltering and multiple hypothesis testing to improve signal
library(ddpca) #Add before cause it requires a version of a package that can be loaded at a lower version in a different package, causing error
library(tidyverse)
#library("arrow")
library(dplyr)
library(readxl)
library(openxlsx)


# COUNT MATRIX: 1. .TSV FILE, 2. ROWNAMES ARE GENE NAMES (SYMBOL) 3. COLNAMES NEED TO COINCIDE WITH A FEATURE FROM METADATA "MATCH_FEATURE"
# FOR UNIT TEST THIS FEATURE IS ANIMAL_ID.
# METADATA: .XLSX DATA FRAME

# PARAMETERS
DirDataRaw = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/raw"
DirData = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/"
DirPipeline = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R"
DirOutput = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/output/231212_Unit_Test"
dir.create(DirOutput)


#PARAM FOR IMPORT_DATA FUNCTION
# Path to Metadata
Metadata_path = paste(DirDataRaw, "Unit_Test_Metadata_V2.xlsx", sep = '/')
# Path to Count matrix
Input_matrix_path = paste(DirDataRaw, "Unit_Test_Count.tsv", sep = '/')
# Path to output file to store results
Output_file_path = DirOutput
# Feature from metadata that matches Count matrix column names
Match_feature = "Animal_id"
# Sample name path to customize sample name
SampleNamepath = paste(DirData,"SampleName.csv", sep = "/")

#PARAM FOR QC_PRENORMALIZATION

# Number of IQR below first quartile to consider an outlier. Default is 2
nMust = 2

#PARAM FOR PRE-FILTERING
# Prevalence fraction filter. A number between 0 and 1.
Prev_perc = 0.1
# PREV uses only Prevalence for pre-filtering, MNB uses prevances AND mixture of negative binomial as pre-filtering, leaving genes that are a mixture of two or more
# negative binomial functions
PRE_FILTER = "PREV"

#PARAM FOR DESEQ_NORM
# Whether normalization is being performed before or after QC_POSTNORMALIZATION. "POST_QCNORM": After QC_POSTNORMALIZATION, "PRE_QCNORM": Before QC_POSTNORMALIZATION
QCNORM = "PRE_QCNORM"
# Path to samples that should be excluded in a customize manner. Add "" if no list wants to be submitted
# See example in "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/outliers_unitTest.csv"
outliers_path = ""
# The Metadata features used to build the Coarse condition to use in DESeq's design matrix. Add Treatment as the first feature.
CoarseConditions = c("Treatment", "Rna_collection_time_hrs")
# Normalization method. Only has 'Standard' otherwise it will stop
METHOD_NORM = 'Standard'
# Use VST transformation after normalization. Either VST_ON or VST_OFF. Otherwise it will stop.
VST_FILTER = "VST_ON"
# Use SVD truncation to get rid of noise from data. 'SVD_ON' otherwise it will continue without implementing filter
SVD_FILTER = 'SVD_OFF'
# Plot PCA per molecule. 'PCA_PLOT' will plot PCA, otherwise it won't
PlotPCA = 'PCA_PLOT'

# PARAM QC-POSTNORMALIZATION
# Outlier filter method. Either 'GENTLE' for mild outlier detection, better for data sets with low sample size per condition, 'STRONG'
# less conservative method for outlier detection, 'NONE' no outlier detection method
OUTLIER_FILTER = "GENTLE_REP"
# Path containing features needed to identify replicates
RepFeatures_path = paste(DirData,"RepFeatures.csv", sep = "/")


#PARAM FOR DEG_FUNCTION
# List of pairwise conditions to be compared from DESeq object. Should be an .xlsx file with two columns named as treat and untreat.
# the elements in the data frame should have the format CoarseConditions[1]_CoarseConditions[2], using the CoarseConditions structure decided in the
# normalization step. The function will only find DEG for molecules present in Metadata file.

List_contrasts_Path = paste(DirData,"List_Contrasts.xlsx", sep = "/")
# DEG methods: Either 'DESeq' or 'T-Test', otherwise it will stop
DEG_Method = 'DESeq'
# Multiple hypothesis testing method. Either BH or High criticism, otherwise it will stop
MH_Method = 'BH'
# Parameter use in high criticism method only. Significance threshold. A number between 0 and 1
AlphaHC = 0.1
# P adjusted value threshold for significance. A number between 0 and 1
padjval = 0.2
# LogFold threshold to use in Volcano plot. A number between 0 and 10000.
LogFoldThrs_VolPlot = 1

# PARAM FOR PATHWAY ENRICHMENT
# Whether the reverse option should be run. Etiher TRUE or FALSE
WITH_REVERSE = TRUE
# Name of the healthy control taken as reference. The DEG between stimulated and this healthy control will be used as a reference to select
# DEG genes from other treatments going in the healthy direction.
REVERSE_DRUG = 'Sham'
# Mostly of in vivo cases there will be only one approximate healthy direction. Select the time to identify this condition
REVERSE_TIME = "72"
# P adjusted value threshold for significance of DEG. A number between 0 and 1
padjval = 0.2
# Log fold change threshold for genes included in pathway enrichment. A value between 0 and 10
LogFoldThrs = 1
# P adjusted value threshold for significance of PATHWAY ENRICHMENT. A number between 0 and 1
Pway_qvalThrs = 0.2
# Organism
ORGANISM = 'Mouse'
# The list of Target candidates from Chemoproteomics
Target_list_path = paste(DirData,"Target_list_Mouse_230922.RData", sep = "/")



# IMPORT DATA -----------------------------------------------------------------------------------
source(paste(DirPipeline, "IMPORT_DATA.R", sep = '/'))
IMPORT_DATA(Metadata_path, Count_path, Match_feature, Output_file_path, SampleNamepath)


# QC-PRENORMALIZATION -----------------------------------------------------------------------------
source(paste(DirPipeline, "QC_PRENORMALIZATION.R", sep = '/'))
QC_PRENORMALIZATION(Output_file_path,nMust)

# PRE-FILTERING -----------------------------------------------------------------------------
source(paste(DirPipeline, "PRE_FILTERING.R", sep = '/'))
PRE_FILTERING(Output_file_path, Prev_perc, PRE_FILTER)

# DESEQ NORMALIZATION -----------------------------------------------------------------------------
source(paste(DirPipeline, "DESEQ_NORM.R", sep = '/'))
DESEQ_NORM(Output_file_path, QCNORM = "PRE_QCNORM", outliers_path = "", CoarseConditions, METHOD_NORM = 'Standard', VST_FILTER = "VST_ON",  SVD_FILTER = 'SVD_OFF', PlotPCA = 'PCA_PLOT',Control_Neg_PCA = "Vehicle", Control_Pos_PCA = "Sham")

# QC-POSTNORMALIZATION-----------------------------------------------------------------------------
source(paste(DirPipeline, "QC_POSTNORMALIZATION.R", sep = '/'))
QC_POSTNORMALIZATION(Output_file_path, OUTLIER_FILTER = OUTLIER_FILTER, RepFeatures_path = RepFeatures_path)

# DESEQ NORMALIZATION -----------------------------------------------------------------------------
source(paste(DirPipeline, "DESEQ_NORM.R", sep = '/'))
CoarseConditions = c("Treatment", "Rna_collection_time_hrs")
outliers_path = ""
DESEQ_NORM(Output_file_path, QCNORM = "POST_QCNORM", outliers_path= outliers_path, CoarseConditions, METHOD_NORM = 'Standard', VST_FILTER = "VST_ON",  SVD_FILTER = 'SVD_OFF', PlotPCA = 'PCA_PLOT', Control_Neg_PCA = "Vehicle", Control_Pos_PCA = "Sham")

# DESEQ DEG -----------------------------------------------------------------------------
source(paste(DirPipeline, "DEG_FUNCTION.R", sep = '/'))
DEG_FUNCTION(Output_file_path, List_contrasts_Path, DEG_Method ,MH_Method, AlphaHC,  padjval  ,LogFoldThrs_VolPlot )

# PATHWAY ENRICHMENT -----------------------------------------------------------------------------
source(paste(DirPipeline, "PWAY_ENRICHMENT.R", sep = '/'))
PWAY_ENRICHMENT(Output_file_path, WITH_REVERSE = TRUE, REVERSE_DRUG = 'Sham', REVERSE_TIME = "72", padjval = 0.2, LogFoldThrs = 1, Pway_qvalThrs = 0.2, ORGANISM = 'Mouse', Target_list_path)







