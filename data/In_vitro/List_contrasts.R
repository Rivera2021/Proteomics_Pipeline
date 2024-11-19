library("readxl")
library(dplyr)
#PARAM FOR IMPORT_DATA FUNCTION
# Path to Metadata

load("~/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0142/Results_Test_Pipeline_241115/DESEQ_NORM/DESeq_Norm.RData")

cellline = 'K562.1'
Metadata =  Metadata %>% filter(Cell_line == cellline)

Mols = unique(Metadata$Treatment)
Mols = Mols[!Mols %in% c('NONE', 'DMSO')]
Mols = "rapamycin"
Df = data.frame(treat = character(), untreat = character())

for(mol in Mols){

    Metadata_mol = Metadata %>% filter(Treatment == mol)

    doses = unique(Metadata_mol$Treatment_conc)
    doses = doses[1:3]
    for(dose in doses){

        Metadata_mol_dose = Metadata_mol %>% filter(Treatment_conc == dose)
        Times = unique(Metadata_mol_dose$Treatment_time_hrs)
        for(Time in Times){

            Metadata_mol_dose_time = Metadata_mol_dose %>% filter(Treatment_time_hrs == Time)
            stims = unique(Metadata_mol_dose_time$Stimulant_used)

            for(stim in stims){

                Metadata_mol_dose_time_stim = Metadata_mol_dose_time %>% filter(Stimulant_used == stim)
                percs = unique(Metadata_mol_dose_time_stim$Reagent_Percentage)

                for(perc in percs){

                    Metadata_mol_dose_time_stim_perc = Metadata_mol_dose_time_stim %>% filter(Reagent_Percentage == perc)
                    Df = Df %>% add_row('treat' = paste(mol,dose,Time, stim,perc,sep = '_'), 'untreat' = paste('DMSO','0',Time, stim,perc, sep = '_'))




                }


            }


        }

    }

}

write.xlsx(Df, file = "~/BULK-TRANSCRIPTOMICS/Cell_line_Experiments/A_2024_0139_AND_A_2024_0140_A_2024_0142_A_2024_0143/Data/List_contrasts_K562.xlsx")
