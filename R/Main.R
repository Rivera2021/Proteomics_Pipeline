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
#PARAM FOR IMPORT_DATA FUNCTION
# Path to Metadata
Metadata_path = '~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/raw/Unit_Test_Metadata_V1.xlsx'
# Path to Count matrix
Input_matrix_path = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/raw/Unit_Test_Count.tsv"
# Path to output file to store results
Output_file_path = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/output/231030_Unit_Test"
# Feature from metadata that matches Count matrix column names
Match_feature = "Animal_id"
# Suffix to delete to shortened name of samples
Suffix_del_name = "Mouse_TNBS_"

#PARAM FOR QC_PRENORMALIZATION
# Number of standard deviations below to consider an optional outlier. Should be higher than 1.
nOpt = 2
# Number of standard deviations below to consider an outlier. Should be higher than 1.
nMust = 4

#PARAM FOR PRE-FILTERING
# Prevalence fraction filter. A number between 0 and 1.
Prev_perc = 0.1
# PREV uses only Prevalence for pre-filtering, MNB uses prevances AND mixture of negative binomial as pre-filtering, leaving genes that are a mixture of two or more
# negative binomial functions
PRE_FILTER = "PREV"

#PARAM FOR DESEQ_NORM
# If a list of samples should be excluded for a certain reason not detected in QC steps
outliers_path = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/outliers_unitTest.csv"
# The Metadata features used to build the Coarse condition to use in DESeq's design matrix. Add Treatment as the first feature.
CoarseConditions = c("Treatment", "Rna_collection_time(hrs)")
# Normalization method. Only has 'Standard' otherwise it will stop
METHOD_NORM = 'Standard'
# Use VST transformation after normalization. Either VST_ON or VST_OFF. Otherwise it will stop.
VST_FILTER = "VST_ON"
# Use SVD truncation to get rid of noise from data. 'SVD_ON' otherwise it will continue without implementing filter
SVD_FILTER = 'SVD_OFF'
# Plot PCA per molecule. 'PCA_PLOT' will plot PCA, otherwise it won't
PlotPCA = 'PCA_PLOT'


#PARAM FOR DEG_FUNCTION
# List of pairwise conditions to be compared from DESeq object. Should be an .xlsx file with two columns named as treat and untreat.
# the elements in the data frame should have the format CoarseConditions[1]_CoarseConditions[2], using the CoarseConditions structure decided in the
# normalization step. The function will only find DEG for molecules present in Metadata file.

List_contrasts_Path = "~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/data/Data_for_pipeline/List_Contrasts.xlsx"

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


# IMPORT DATA -----------------------------------------------------------------------------------
source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/IMPORT_DATA.R")
IMPORT_DATA(Metadata_path, Count_path, Match_feature, Output_file_path, Suffix_del_name)


# QC-PRENORMALIZATION -----------------------------------------------------------------------------
source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/QC_PRENORMALIZATION.R")
QC_PRENORMALIZATION(Output_file_path, nOpt, nMust)

# PRE-FILTERING -----------------------------------------------------------------------------
source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/PRE_FILTERING.R")
PRE_FILTERING(Output_file_path, Prev_perc, PRE_FILTER)

# DESEQ PRENORMALIZATION -----------------------------------------------------------------------------
source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/DESEQ_NORM.R")
DESEQ_NORM(Output_file_path, outliers_path, CoarseConditions, METHOD_NORM = 'Standard', VST_FILTER = "VST_ON",  SVD_FILTER = 'SVD_OFF', PlotPCA = 'PCA_PLOT')
# DESEQ DEG -----------------------------------------------------------------------------
source("~/BULK-TRANSCRIPTOMICS/Transcriptomics_Pipeline/R/DEG_FUNCTION.R")
DEG_FUNCTION(Output_file_path, List_contrasts_Path, DEG_Method = 'DESeq',MH_Method = 'BH', AlphaHC = 0.1,  padjval = 0.2 ,LogFoldThrs_VolPlot = 1)








