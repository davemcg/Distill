######################
# merge models from build_*PaC.R
#######################
library(randomForest)
library(tidyverse)

model_merger <- function(pattern='VPaC__6mtry_v2*', mtry=6, pac = 'VPaC_pieces', type='VPaC', version='v1'){
  files <- list.files(paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/', pac), 
                      pattern = pattern, 
                      full.names = TRUE)
  for (i in files){load(i)}
  ls()[grepl('VPaC', ls())]
  
  pac_list <- list()
  for (i in ls()[grepl('VPaC', ls())]){pac_list[[i]] <- get(i)}
  PaC <- do.call(randomForest::combine, pac_list)
  
  
  name = paste0(type, '_', mtry, 'mtry', '_', version)
  assign(name, PaC)
  ##############################
  # SAVE MODEL
  ###############################
  save(list=name, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/', name, '.Rdata'))
}
# 
# #VPaC
# model_merger(pattern='VPaC__9mtry_v3*', mtry=9, pac='VPaC_pieces', type='VPaC')
# v

#VPaC
model_merger(pattern='VPaC__06mtry_v7_*', mtry=6, pac='VPaC_pieces', type='VPaC', version = 'v8')
model_merger(pattern='VPaC__09mtry_v7_*', mtry=9, pac='VPaC_pieces', type='VPaC', version = 'v8')
model_merger(pattern='VPaC__12mtry_v7_*', mtry=12, pac='VPaC_pieces', type='VPaC', version = 'v8')


