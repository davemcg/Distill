
library(tidyverse)
library(data.table)
library(dummies)
library(caret)
library(mlbench)
library(parallel)
library(doParallel)
library(MLmetrics)

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
most_imp_predictors <- c('is_lof','impact_severity','`DiseaseClass_-1`','mis_z','ccr_pct_v1','cadd_phred','phylop_100way','n_mis','revel','fitcons_float','precessive','n_lof','`DiseaseClass_Stargardt,RD`','m_cap_rankscore','dann','DiseaseClass_RD','vest3_rankscore','n_syn','pnull','pli','lof_z','fathmm_mkl_coding_rankscore','an_exac_all','eigen_pc_raw_rankscore','gerp_elements','mutationassessor_score_rankscore','mpc','metasvm_rankscore','polyphen_score','metalr_rankscore','lrt_converted_rankscore','genocanyon_score_rankscore','mutationtaster_converted_rankscore','gno_an_popmax','grantham','max_aaf_all','ac_exac_all','fathmm_converted_rankscore','aaf_esp_all','sift_score','ac_exac_sas', 'linsight')
most_imp_predictors_no_disease_class <- c('is_lof','impact_severity','mis_z','ccr_pct_v1','cadd_phred','phylop_100way','n_mis','revel','fitcons_float','precessive','n_lof','m_cap_rankscore','dann','vest3_rankscore','n_syn','pnull','pli','lof_z','fathmm_mkl_coding_rankscore','an_exac_all','eigen_pc_raw_rankscore','gerp_elements','mutationassessor_score_rankscore','mpc','metasvm_rankscore','polyphen_score','metalr_rankscore','lrt_converted_rankscore','genocanyon_score_rankscore','mutationtaster_converted_rankscore','gno_an_popmax','grantham','max_aaf_all','ac_exac_all','fathmm_converted_rankscore','aaf_esp_all','sift_score','ac_exac_sas','linsight')

###########################################
# multi processing
##########################################
cluster <- makeCluster(10) 
registerDoParallel(cluster)

##############################################
# BUILD MODELS!!!!!!!!!!!
#############################################


rfFit_OVPaC <- caret::train(Status ~ ., data=model_data$ML_set__eye_dummy_TT$train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                            method = "rf", metric='F', ntree=1501,
                            trControl = fitControl_min)


##############################
# SAVE MODEL
###############################

save(rfFit_OVPaC, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/rfFit_OVPaC.Rdata')
