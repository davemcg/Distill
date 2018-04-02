
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
# rfFit_all <- caret::train(Status ~ ., data=train_set %>% select(-pos_id, -Source), 
#                       method = "rf", metric='F',
#                       trControl = fitControl_min)
# # use the first rf model to pick the semi-useful predictors and limit the models to these
# most_imp_predictors <- varImp(rfFit_all)$importance  %>% rownames_to_column('Predictors') %>% arrange(-Overall) %>% filter(Overall > 4) %>% pull(Predictors)
# # variant with no disease class predictors
# most_imp_predictors_no_disease_class <- most_imp_predictors[!grepl('DiseaseClass', most_imp_predictors)]
# 
# model_run$rfFit_all <- rfFit_all
# model_run$most_imp_predictors <- most_imp_predictors
# model_run$most_imp_predictors_no_disease_class <- most_imp_predictors_no_disease_class
# save(model_run, file='model_run__2018_03_29.Rdata')

rfFit_OVPaC <- caret::train(Status ~ ., data=model_data$ML_set__eye_dummy_TT$train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
                      method = "rf", metric='F', ntree=1501,
                      trControl = fitControl_min)

# rfFit_VPaC <- caret::train(Status ~ ., data=model_data$ML_set__general_dummy_TT$train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                             method = "rf", metric='F', ntree=10,
#                             trControl = fitControl_min)

# bglmFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                         method = "bayesglm", metric='F',
#                         trControl = fitControl_min)
# 
# bglmFit_noDC <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
#                              method = "bayesglm", metric='F',
#                              trControl = fitControl_min)

# glmFit <- caret::train(Status ~ ., data=validate_set %>% select_(.dots=c('Status','pathUK','pathC')), 
#                        method = "glm", metric='F',
#                        trControl = fitControl_min)
# 
# glmFit_noDC <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
#                             method = "glm", metric='F',
#                             trControl = fitControl_min)

# glmboostFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                             method = "glmboost", metric='F',
#                             trControl = fitControl_min,
#                             preProcess = c('center','scale'))

# LogitBoostFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                               method = "LogitBoost", metric='F',
#                               trControl = fitControl_min,
#                               preProcess = c('center','scale'))

# avNNetFit <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors)), 
#                           method = "avNNet", metric='F',
#                           trControl = fitControl_min,
#                           preProcess = c('center','scale'))
# 
# avNNetFit_noDC <- caret::train(Status ~ ., data=train_set %>% select_(.dots=c('Status',most_imp_predictors_no_disease_class)), 
#                                method = "avNNet", metric='F',
#                                trControl = fitControl_min,
#                                preProcess = c('center','scale'))

# tossed, terrible performance
# xgbTreeFit <- caret::train(Status ~ ., data=train_set %>% select(-variant_id, -Source), 
#                      
#                       method = "xgbTree",  metric='AUC',
#                       trControl = fitControl)

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

##############################
# SAVE MODELS
###############################

for (i in ls()[grepl('Fit',ls())]) {model_run[[i]] <- get(i)}
save(model_run, file='model_run__2018_03_28.Rdata')
