#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

model <- args[1]
cores <- as.numeric(args[2])

library(tidyverse)
library(data.table)
#library(dummies)
library(caret)
library(mlbench)
library(parallel)
library(doParallel)
library(ModelMetrics)
library(PRROC)

#load('/Volumes/Arges/PROJECTS/mcgaughey/eye_var_Pathogenicity/clean_data/model_data_2018_07_13.Rdata')
#load('/Volumes/Arges/PROJECTS/mcgaughey/eye_var_Pathogenicity/clean_data/assess_2018_07_17.Rdata')

########################################
# Set up fitcontrols for caret modeling
#########################################
# CV on rf seems to overfit
fitControl_min <- trainControl(
  classProbs=T,
  savePredictions = T,
  allowParallel = T,
  summaryFunction = prSummary,
  returnData = T,
  preProcOptions = c('center','scale'))

# fitControl_min_CS <- trainControl(
#   classProbs=T,
#   savePredictions = T,
#   allowParallel = T,
#   summaryFunction = prSummary,
#   returnData = T,
#   preProcOptions = c('center','scale'))

# fitControl <- trainControl(## 5-fold CV
#   method = "cv",
#   number = 5,
#   classProbs=T,
#   savePredictions = T,
#   allowParallel = T,
#   summaryFunction = prSummary,
#   returnData = T)
######################################
# Predictors
#######################################
most_imp_predictors <- c('ccr_pct_v1','cadd_raw','vest3_rankscore','cadd_phred','mis_z','pli','lof_z','phylop_100way','revel','hapmap2','hapmap1','n_mis','epilogos_quies','n_lof','precessive','pnull','adj_exp_lof','adj_exp_syn','dann','adj_exp_mis','syn_z','n_syn','epilogos_txwk','fitcons','m_cap_score','m_cap_rankscore','eigen_phred','eigen_raw','epilogos_tx','eigen_pc_raw_rankscore','epilogos_reprpcwk','fathmm_mkl_coding_rankscore','metalr_score','fathmm_mkl_coding_score','metalr_rankscore','impact_severity','metasvm_rankscore','metasvm_score','epilogos_enh','genocanyon_score','sigmaaf_lof_0001','sigmaaf_lof_01','sigmaaf_missense_0001','sigmaaf_missense_01')

most_imp_predictors_expand <- c('ccr_pct_v1','cadd_raw','vest3_rankscore','cadd_phred','mis_z','pli','lof_z','phylop_100way','revel','hapmap2','hapmap1','n_mis','epilogos_quies','n_lof','precessive','pnull','adj_exp_lof','adj_exp_syn','dann','adj_exp_mis','syn_z','n_syn','epilogos_txwk','fitcons','m_cap_score','m_cap_rankscore','eigen_phred','eigen_raw','epilogos_tx','is_lof','eigen_pc_raw_rankscore','epilogos_reprpcwk','fathmm_mkl_coding_rankscore','metalr_score','fathmm_mkl_coding_score','metalr_rankscore','impact_severity','metasvm_rankscore','metasvm_score','epilogos_enh','genocanyon_score','fathmm_converted_rankscore','mpc','epilogos_enhg','af_exac_all','epilogos_reprpc','max_aaf_all','mutationassessor_score','gerp','polyphen_score','gerp_elements','mutationassessor_score_rankscore','stam_mean','an_exac_all','af_exac_nfe','provean_converted_rankscore','an_exac_nfe','lrt_score','lrt_omega','grantham','lrt_converted_rankscore','genocanyon_score_rankscore','an_exac_afr','an_exac_amr','an_exac_sas','epilogos_het','ac_exac_all','linsight','gno_an_popmax','exac_num_het','an_exac_eas','gno_an_all','ac_exac_nfe','mutationtaster_converted_rankscore','an_exac_oth','an_exac_fin','gno_an_nfe','gno_af_all','gno_an_afr','epilogos_tssaflnk','gno_af_popmax','epilogos_znf','segway_sum_score','aaf_esp_ea','epilogos_txflnk','provean_score','segway_mean_score','epilogos_tss','aaf_esp_all','af_exac_amr','gno_af_nfe','epilogos_enhbiv','af_exac_sas','sift_score','fathmm_score','ac_exac_amr','aaf_esp_aa','gno_ac_all','gno_af_afr','ac_exac_sas','af_exac_eas','gno_an_fin','af_exac_afr','gno_an_eas','gno_an_oth','gno_ac_nfe','gno_ac_popmax','ac_exac_eas','ac_exac_afr','epilogos_tssbiv','gno_ac_afr','vest3_score','sigmaaf_lof_0001','sigmaaf_lof_01','sigmaaf_missense_0001','sigmaaf_missense_01')

###########################################
# multi processing
##########################################
print(paste(cores, "cores loading"))
cluster <- makeCluster(cores) 
registerDoParallel(cluster)

#############################################
# cut data down to speed up evaluation
##############################################
# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_data_2018_07_13.Rdata')
# load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/assess_2018_07_17.Rdata')
# set.seed(52349)
# 
# train_Path <- model_data$ML_set__general_TT$train_set %>% filter(Status == 'Pathogenic') %>% sample_frac(0.2)
# train_NotPath <- model_data$ML_set__general_TT$train_set %>% filter(Status == 'NotPathogenic') %>% sample_n(nrow(train_Path) * 25)
# train_data <- bind_rows(train_Path, train_NotPath)
# train_data <- train_data %>% sample_n(nrow(train_data))
# test_data <- assess_set %>% filter(DataSet == 'SuperGrimm') %>% sample_frac(0.1)
# model_toy_data <- list()
# model_toy_data$train_data <- train_data
# model_toy_data$test_data <- test_data
# save(model_toy_data, file='/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_toy_data_2018_07_18.Rdata')

load('/data/mcgaugheyd/projects/nei/mcgaughey/eye_var_Pathogenicity/clean_data/model_toy_data_2018_07_18.Rdata')
train_data <- model_toy_data$train_data
test_data <- model_toy_data$test_data

##############################################
# BUILD MODELS!!!!!!!!!!!
#############################################
print(paste("Starting",model,"build"))
model <- caret::train(Status ~ ., data = train_data %>% select(one_of(c('Status', most_imp_predictors))),
                      method = model, metric='F',
                      trControl = fitControl_min)

# glmFit <- caret::train(Status ~ ., data = train_data %>% select(one_of(c('Status', most_imp_predictors))),
#                        method = "glm", metric='F',
#                        trControl = fitControl_min)
# 
# rfFit <- caret::train(Status ~ ., data = train_data %>% select(one_of(c('Status',most_imp_predictors_expand))), 
#                       method = "rf", metric='F', ntree=50,
#                       trControl = fitControl_min)
# 
# xgbTreeFit <- caret::train(Status ~ ., data=train_data %>% select(one_of(c('Status', most_imp_predictors))),
#                            method = "xgbTree",  metric='F',
#                            trControl = fitControl_min)
# # 
# # xgbDARTFit <- caret::train(Status ~ ., data=train_data %>% select(one_of(c('Status', most_imp_predictors))),
# #                            method = "xgbDART",  metric='F',
# #                            trControl = fitControl_min)
# 
# xgbLinearFit <- caret::train(Status ~ ., data=train_data %>% select(one_of(c('Status', most_imp_predictors))),
#                              method = "xgbLinear",  metric='F',
#                              trControl = fitControl_min)
# 
# lssvmPolyFit <- caret::train(Status ~ ., data=train_data %>% select(one_of(c('Status', most_imp_predictors))),
#                              method = "lssvmPoly",  metric='F',
#                              trControl = fitControl_min)
# 
# svmPolyFit <- caret::train(Status ~ ., data=train_data %>% select(one_of(c('Status', most_imp_predictors))),
#                            method = "svmPoly",  metric='F',
#                            trControl = fitControl_min)
# 
# # adaboostFit <- caret::train(Status ~ ., data=train_data %>% select(one_of(c('Status', most_imp_predictors))),
# #                             method = "adaboost",  metric='F',
# #                             trControl = fitControl_min)
# 
# roccFit <-  caret::train(Status ~ ., data=train_data %>% select(one_of(c('Status', most_imp_predictors))),
#                          method = "rocc",  metric='F',
#                          trControl = fitControl_min)
# 
# xyfFit <- caret::train(Status ~ ., data=train_data %>% select(one_of(c('Status', most_imp_predictors))),
#                        method = "xyf",  metric='F',
#                        trControl = fitControl_min)


# stepLDAFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                            method = "stepLDA",  metric='AUC',
#                            trControl = fitControl_min_CS)

# tossed, terrible performance
# naive_bayesFit <- caret::train(Status ~ ., data=train_set %>% select(-variant_id, -Source), 
#                      
#                       method = "naive_bayes",  metric='AUC',
#                       trControl = fitControl)

# dnnFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                        method = "dnn",  metric='AUC',
#                        trControl = fitControl_min,
#                        preProcess = c('center','scale'))
# 
# dnnFit_noDC <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
#                        method = "dnn",  metric='AUC',
#                        trControl = fitControl_min,
#                        preProcess = c('center','scale'))


# monmlpFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                           method = "monmlp",  metric='AUC',
#                           trControl = fitControl_min)


# monmlpFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                           method = "monmlp",  metric='F',
#                           trControl = fitControl_min)

# mlpKerasDropoutFit <- caret::train(x = train_set %>% select(-variant_id, -Source, -Status),
#                                  y= train_set$Status,
#                      
#                       method = "mlpKerasDropout",  metric='F',
#                       trControl = fitControl)

# svmLinearWeightsFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                                     method = "svmLinearWeights",  metric='Precision',
#                                     trControl = fitControl_min,
#                                     preProcess = c('center','scale'))

# caddFit <- caret::train(Status ~ ., data=train_set %>% select(cadd_phred, Status), 
#                         method = "glm",  metric='AUC',
#                         trControl = fitControl_min)
# 
# cadd_DC_Fit <- caret::train(Status ~ ., data=train_set %>% select(cadd_phred, Status, `DiseaseClass_-1`), 
#                             method = "glm",  metric='AUC',
#                             trControl = fitControl_min)

# revelFit <- caret::train(Status ~ ., data=train_set %>% select(revel, Status), 
#                          method = "glm", metric='AUC',
#                          trControl = fitControl_min)
# 
# revel_DC_Fit <- caret::train(Status ~ ., data=train_set %>% select(revel, Status, `DiseaseClass_-1`), 
#                              method = "glm", metric='AUC',
#                              trControl = fitControl_min)

# dannFit <- caret::train(Status ~ ., data=train_set %>% select(dann, Status), 
#                         method = "glm",
#                         trControl = fitControl_min)
# 
# dann_DC_Fit <- caret::train(Status ~ ., data=train_set %>% select(dann, Status, `DiseaseClass_-1`), 
#                             method = "glm",
#                             trControl = fitControl_min)

################################
# Quick Assess
################################
cm_maker <- function(predictor = 'cadd_phred', data, cutoff=0.5, mode = 'prec_recall') {
  if (class(predictor)!='character'){
    print("Running in predictor is a model mode")
    new_predictions <- predict(predictor, data, type='prob') %>% data.frame() %>% 
      mutate(Answers = data$Status, Prediction = case_when(Pathogenic > cutoff ~ 'Pathogenic', TRUE ~ 'NotPathogenic'))
    new_predictions <- new_predictions %>% mutate(preds = case_when(Prediction == 'Pathogenic' ~ 1,
                                                                    TRUE ~ 0),
                                                  actuals = case_when(Answers == 'Pathogenic' ~ 1,
                                                                      TRUE ~ 0))
    out <- caret::confusionMatrix(data = as.factor(new_predictions$Prediction), reference = as.factor(new_predictions$Answers), mode= mode)
    out$MCC <- mcc(new_predictions$preds, new_predictions$actuals, cutoff=cutoff)
  } else {
    print("Running in predictor is a precomputed column in data mode")
    new_predictions <- data 
    new_predictions$Prediction <- 'NotPathogenic'
    new_predictions[(new_predictions[,predictor] > cutoff), 'Prediction'] <- "Pathogenic"
    new_predictions <- new_predictions %>% mutate(preds = case_when(Prediction == 'Pathogenic' ~ 1,
                                                                    TRUE ~ 0),
                                                  actuals = case_when(Status == 'Pathogenic' ~ 1,
                                                                      TRUE ~ 0))
    out <- caret::confusionMatrix(data = as.factor(new_predictions$Prediction), reference = as.factor(new_predictions$Status), mode= mode)
    out$MCC <- mcc(new_predictions$preds, new_predictions$actuals, cutoff=cutoff)
  }
  out
}

print(model$method)
metrics <- cm_maker(model, test_data)
print(metrics$MCC)
print(metrics$overall)

##############################
# SAVE MODELS
###############################

# for (i in ls()[grepl('Fit',ls())]) {model_run[[i]] <- get(i)}
# save(model_run, file='model_run__2018_03_28.Rdata')
