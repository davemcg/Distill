######################
# merge models from build_*PaC.R
#######################
library(randomForest)
library(tidyverse)

#VPaC
vpac_files <- list.files('/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_pieces/', pattern = '*Rdata', full.names = TRUE)
for (i in vpac_files){load(i)}
ls()[grepl('VPaC', ls())]

vpac_list <- list()
for (i in ls()[grepl('VPaC', ls())]){vpac_list[[i]] <- get(i)}
VPaC <- do.call(randomForest::combine, vpac_list)
save(VPaC, file ='/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/VPaC_RF.Rdata')

#OVPaC
ovpac_files <- list.files('/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/OVPaC_pieces/', pattern = '*Rdata', full.names = TRUE)
for (i in ovpac_files){load(i)}
ls()[grepl('OVPaC', ls())]

ovpac_list <- list()
for (i in ls()[grepl('OVPaC', ls())]){ovpac_list[[i]] <- get(i)}
OVPaC <- do.call(randomForest::combine, ovpac_list)
save(OVPaC, file ='/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/OVPaC_RF.Rdata')

