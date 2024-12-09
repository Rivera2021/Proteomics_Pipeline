# GENERATE THE CONTRAST LIST REQUIRED BY THE PIPELINE

library("openxlsx")
library("readxl")
library(dplyr)
#PARAM FOR IMPORT_DATA FUNCTION
# Path to Metadata

load("~/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/EXP000072/Results_241102/DESEQ_NORM/DESeq_Norm.RData")
#CoarseCondition:  ["Cell_line","Treatment","Treatment_conc", "Treatment_time_hrs", "Stimulant_used"]
cellline = 'K562.3'
Metadata =  Metadata %>% dplyr::filter(Cell_line == cellline)

Mols = unique(Metadata$Treatment)
Mols = Mols[!Mols %in% c('NONE', 'DMSO')]

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
                Df = Df %>% dplyr::add_row('treat' = paste(cellline, mol,dose,Time, stim,sep = '_'), 'untreat' = paste(cellline,'DMSO','0',Time, stim, sep = '_'))



            }


        }

    }

}

write.xlsx(Df, file = "~/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/EXP000072/Data/List_contrasts_K562.xlsx")

