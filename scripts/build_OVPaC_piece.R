########################
# Build VPaC
# generic pathogenicity score
#########################


############################
# Only makes 33 trees
# Run 50 times and combine the 50 33-tree forests with `randomForest::combine`
#############################

############################
# give 8GB for biowulf2 run
############################

library(tidyverse)
library(randomForest)
library(caret)

load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_data.Rdata')

########################################
# Set up fitcontrols for caret modeling
#########################################
# CV on rf seems to overfit
fitControl_min <- trainControl(
  classProbs=T,
  savePredictions = T,
  allowParallel = T,
  summaryFunction = prSummary,
  returnData = T)

######################################
# Predictors
#######################################
most_imp_predictors <- c('is_lof','impact_severity','DiseaseClass','mis_z','ccr_pct_v1','cadd_phred','phylop_100way','n_mis','revel','fitcons_float','precessive','n_lof','m_cap_rankscore','dann','vest3_rankscore','n_syn','pnull','pli','lof_z','fathmm_mkl_coding_rankscore','an_exac_all','eigen_pc_raw_rankscore','gerp_elements','mutationassessor_score_rankscore','mpc','metasvm_rankscore','polyphen_score','metalr_rankscore','lrt_converted_rankscore','genocanyon_score_rankscore','mutationtaster_converted_rankscore','gno_an_popmax','grantham','max_aaf_all','ac_exac_all','fathmm_converted_rankscore','aaf_esp_all','sift_score','ac_exac_sas', 'linsight')
most_imp_predictors_no_disease_class <- c('is_lof','impact_severity','mis_z','ccr_pct_v1','cadd_phred','phylop_100way','n_mis','revel','fitcons_float','precessive','n_lof','m_cap_rankscore','dann','vest3_rankscore','n_syn','pnull','pli','lof_z','fathmm_mkl_coding_rankscore','an_exac_all','eigen_pc_raw_rankscore','gerp_elements','mutationassessor_score_rankscore','mpc','metasvm_rankscore','polyphen_score','metalr_rankscore','lrt_converted_rankscore','genocanyon_score_rankscore','mutationtaster_converted_rankscore','gno_an_popmax','grantham','max_aaf_all','ac_exac_all','fathmm_converted_rankscore','aaf_esp_all','sift_score','ac_exac_sas','linsight')


##############################################
# BUILD MODEL!!!!!!!!!!!
#############################################

rand_num <- as.integer(paste(sample(0:9, 6, replace=F), collapse = ''))
set.seed(rand_num)

rfFit_OVPaC <- randomForest(Status ~ ., data=model_data$ML_set__eye_TT$train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                            ntree=20,
                            mtry=10,
                            importance = TRUE,
                            norm.votes = FALSE)

rand_name1 = paste0('OVPaC__10mtry_v2', rand_num)
assign(rand_name1, rfFit_OVPaC)


rfFit_OVPaC <- randomForest(Status ~ ., data=model_data$ML_set__eye_TT$train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                            ntree=20,
                            mtry=6,
                            importance = TRUE,
                            norm.votes = FALSE)

rand_name2 = paste0('OVPaC__6mtry_v2', rand_num)
assign(rand_name2, rfFit_OVPaC)


rfFit_OVPaC <- randomForest(Status ~ ., data=model_data$ML_set__eye_TT$train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                            ntree=20,
                            mtry=3,
                            importance = TRUE,
                            norm.votes = FALSE)

rand_name3 = paste0('OVPaC__3mtry_v2', rand_num)
assign(rand_name3, rfFit_OVPaC)


##############################
# SAVE MODEL
###############################
save(list = rand_name1, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/OVPaC_pieces/', rand_name1, '.Rdata'))
save(list = rand_name2, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/OVPaC_pieces/', rand_name2, '.Rdata'))
save(list = rand_name3, file=paste0('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/OVPaC_pieces/', rand_name3, '.Rdata'))
