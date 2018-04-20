######################
# merge models from build_*PaC.R
#######################
library(randomForest)
library(tidyverse)

model_merger <- function(pattern='VPaC__6mtry_v2*', mtry=6, pac = 'VPaC_pieces', type='VPaC'){
  files <- list.files(paste0('/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/', pac), pattern = pattern, full.names = TRUE)
  for (i in files){load(i)}
  ls()[grepl('VPaC', ls())]
  
  pac_list <- list()
  for (i in ls()[grepl('VPaC', ls())]){pac_list[[i]] <- get(i)}
  PaC <- do.call(randomForest::combine, pac_list)
  
  
  name = paste0(type, '_', mtry, 'mtry')
  assign(name, PaC)
  ##############################
  # SAVE MODEL
  ###############################
  save(list=name, file=paste0('/Volumes/data/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/', name, '.Rdata'))
}

#VPaC
model_merger(pattern='VPaC__3mtry_v2*', mtry=3, pac='VPaC_pieces', type='VPaC')
model_merger(pattern='VPaC__6mtry_v2*', mtry=6, pac='VPaC_pieces', type='VPaC')
model_merger(pattern='VPaC__10mtry_v2*', mtry=10, pac='VPaC_pieces', type='VPaC')

#OVPaC
model_merger(pattern='OVPaC__3mtry_v2*', mtry=3, pac='OVPaC_pieces', type='OVPaC')
model_merger(pattern='OVPaC__6mtry_v2*', mtry=6, pac='OVPaC_pieces', type='OVPaC')
model_merger(pattern='OVPaC__10mtry_v2*', mtry=10, pac='OVPaC_pieces', type='OVPaC')


