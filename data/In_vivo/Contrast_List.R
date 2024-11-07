# GENERATE THE CONTRAST LIST REQUIRED BY THE PIPELINE

library("openxlsx")
library("readxl")
library(dplyr)
#PARAM FOR IMPORT_DATA FUNCTION
# Path to Metadata

load("~/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Results_240929/DESEQ_NORM/DESeq_Norm.RData")


Mols = unique(Metadata$Treatment)
Mols = Mols[!Mols %in% c('NONE', 'Vehicle-DMSO')]

Df = data.frame(treat = character(), untreat = character())

for(mol in Mols){

    Metadata_mol = Metadata %>% dplyr::filter(Treatment == mol)

    doses = unique(Metadata_mol$Treatment_conc)
    for(dose in doses){

        Metadata_mol_dose = Metadata_mol %>% dplyr::filter(Treatment_conc == dose)
        Times = unique(Metadata_mol_dose$Treatment_time_hrs)
        for(Time in Times){

            Metadata_mol_dose_time = Metadata_mol_dose %>% dplyr::filter(Treatment_time_hrs == Time)
            stims = unique(Metadata_mol_dose_time$Stimulant_used)

            for(stim in stims){

                Metadata_mol_dose_time_stim = Metadata_mol_dose_time %>% dplyr::filter(Stimulant_used == stim)
                Df = Df %>% dplyr::add_row('treat' = paste(mol,dose,Time,stim,sep = '_'), 'untreat' = paste('Vehicle-DMSO','0',Time,stim, sep = '_'))



            }


        }

    }

}

Df = Df %>% dplyr::add_row('treat' = "NONE_0_72_NoStim"  , 'untreat' = "Vehicle-DMSO_0_6_Stim")
Df = Df %>% dplyr::add_row('treat' = "NONE_0_72_NoStim"  , 'untreat' = "Vehicle-DMSO_0_12_Stim")
Df = Df %>% dplyr::add_row('treat' = "NONE_0_72_NoStim"  , 'untreat' = "Vehicle-DMSO_0_24_Stim")
Df = Df %>% dplyr::add_row('treat' = "NONE_0_72_NoStim"  , 'untreat' = "Vehicle-DMSO_0_72_Stim")

write.xlsx(Df, file = "~/BULK-TRANSCRIPTOMICS/InVivo/EXP000033/Data/List_contrasts_Invivo.xlsx")


